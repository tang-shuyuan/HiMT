"""
view.py - 纯 Python 实现的 Bandage 风格 GFA 图谱可视化渲染器

功能概述：
  1. 解析 GFA（Graphical Fragment Assembly）格式的基因组组装图文件
  2. 使用力导向布局算法（FM³ - Fast Multipole Multilevel Method）自动计算图的节点布局
  3. 将组装图渲染为 PNG 或 SVG 图片

核心流程：
  GFA 文件 → 解析为 DeBruijn 图 → 构建 OGDF 布局图 → 力导向布局计算 → 渲染输出图片

主要模块：
  - 颜色/几何工具：Point2D, RgbaColor, 颜色转换等基础数据结构
  - 图数据结构：DeBruijnNode, DeBruijnEdge, AssemblyGraph（DeBruijn 图模型）
  - GFA 解析器：从 GFA 文件读取节点（S 行）和边（L 行）
  - 布局引擎：GraphLayoutWorker（FM³ 多层级力导向布局算法）
  - 渲染器：render_png / render_svg（将布局结果绘制为图片）
"""

from __future__ import annotations

import argparse
import datetime
import math
import random
import time
from collections import deque
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape


# ==================== 常量定义 ====================

EPSILON = 1e-9      # 浮点数比较的极小阈值，避免除零等问题
TWO_PI = math.pi * 2.0  # 2π，用于角度计算


# ==================== 基础数据结构 ====================

@dataclass(frozen=True)
class Point2D:
    """二维点坐标（不可变）"""
    x: float
    y: float


@dataclass(frozen=True)
class RgbaColor:
    """RGBA 颜色（不可变），取值范围 0-255"""
    r: int
    g: int
    b: int
    a: int = 255  # 透明度，255 表示完全不透明

    def hex(self) -> str:
        """返回十六进制颜色字符串，如 #ff0000"""
        return f"#{self.r:02x}{self.g:02x}{self.b:02x}"

    def opacity(self) -> float:
        """返回透明度比例值（0.0 ~ 1.0）"""
        return self.a / 255.0

    def to_pillow(self) -> tuple[int, int, int, int]:
        """转换为 Pillow 图像库使用的 RGBA 元组"""
        return (self.r, self.g, self.b, self.a)


# 预定义常用颜色
BLACK = RgbaColor(0, 0, 0, 255)
WHITE = RgbaColor(255, 255, 255, 255)
LABEL_HALO_COLOUR = RgbaColor(255, 255, 255, 230)  # 标签光晕色（白色半透明，用于文字描边）


# ==================== 颜色工具函数 ====================

def clamp255(channel: float) -> int:
    """将 [0.0, 1.0] 范围的浮点通道值夹紧并转换为 [0, 255] 的整数"""
    return max(0, min(255, int(round(channel * 255.0))))


def hue_to_rgb(p: float, q: float, t: float) -> float:
    """HSL → RGB 转换的辅助函数，处理单个颜色通道"""
    if t < 0.0:
        t += 1.0
    if t > 1.0:
        t -= 1.0
    if t < 1.0 / 6.0:
        return p + (q - p) * 6.0 * t
    if t < 1.0 / 2.0:
        return q
    if t < 2.0 / 3.0:
        return p + (q - p) * (2.0 / 3.0 - t) * 6.0
    return p


def hsl_to_rgb(hue: float, saturation: float, lightness: float, alpha: int) -> RgbaColor:
    """
    将 HSL 颜色空间转换为 RGBA 颜色
    - hue: 色相 [0.0, 1.0]
    - saturation: 饱和度 [0.0, 1.0]
    - lightness: 亮度 [0.0, 1.0]
    - alpha: 不透明度 [0, 255]
    用于为图中每个节点生成随机颜色
    """
    if saturation == 0.0:
        # 饱和度为 0 时为灰色
        r = lightness
        g = lightness
        b = lightness
    else:
        q = lightness * (1.0 + saturation) if lightness < 0.5 else lightness + saturation - lightness * saturation
        p = 2.0 * lightness - q
        r = hue_to_rgb(p, q, hue + 1.0 / 3.0)
        g = hue_to_rgb(p, q, hue)
        b = hue_to_rgb(p, q, hue - 1.0 / 3.0)
    return RgbaColor(clamp255(r), clamp255(g), clamp255(b), alpha)


def parse_color(value: str) -> RgbaColor:
    """
    解析颜色字符串为 RgbaColor 对象
    支持命名颜色（如 "black", "red"）和十六进制格式（如 "#ff0000", "ff0000", "ffff0000"）
    """
    normalized = value.strip().lower()
    named = {
        "black": BLACK,
        "white": WHITE,
        "red": RgbaColor(255, 0, 0),
        "green": RgbaColor(0, 255, 0),
        "blue": RgbaColor(0, 0, 255),
        "yellow": RgbaColor(255, 255, 0),
        "orange": RgbaColor(255, 200, 0),
        "cyan": RgbaColor(0, 255, 255),
        "magenta": RgbaColor(255, 0, 255),
        "pink": RgbaColor(255, 175, 175),
        "gray": RgbaColor(128, 128, 128),
        "grey": RgbaColor(128, 128, 128),
        "darkgray": RgbaColor(64, 64, 64),
        "darkgrey": RgbaColor(64, 64, 64),
        "lightgray": RgbaColor(192, 192, 192),
        "lightgrey": RgbaColor(192, 192, 192),
    }
    if normalized in named:
        return named[normalized]
    hex_value = normalized[1:] if normalized.startswith("#") else normalized
    if hex_value.startswith("0x"):
        hex_value = hex_value[2:]
    if len(hex_value) == 6:
        # 6 位十六进制：RRGGBB，不透明
        return RgbaColor(int(hex_value[0:2], 16), int(hex_value[2:4], 16), int(hex_value[4:6], 16), 255)
    if len(hex_value) == 8:
        # 8 位十六进制：AARRGGBB，含透明度
        argb = int(hex_value, 16)
        return RgbaColor((argb >> 16) & 0xFF, (argb >> 8) & 0xFF, argb & 0xFF, (argb >> 24) & 0xFF)
    raise ValueError(f"Unsupported font color: {value}. Expected a name like black or a hex color like #000000")


# ==================== 几何/数学工具函数 ====================

def format_number(value: float) -> str:
    """将浮点数格式化为保留 3 位小数的字符串，用于 SVG 输出"""
    return f"{value:.3f}"


def format_trimmed_depth(value: float) -> str:
    """
    格式化深度值用于标签显示
    - 如果接近整数则显示整数（如 5.0 → "5"）
    - 否则保留 1 位小数（如 5.3 → "5.3"）
    """
    if abs(value - round(value)) < 0.05:
        return str(int(round(value)))
    return f"{value:.1f}"


def distance_points(a: Point2D, b: Point2D) -> float:
    """计算两点之间的欧氏距离"""
    return math.hypot(a.x - b.x, a.y - b.y)


def angle_vectors(qx: float, qy: float, rx: float, ry: float) -> float:
    """
    计算两个向量之间的夹角（弧度），返回值范围 [0, 2π)
    - (qx, qy): 第一个向量
    - (rx, ry): 第二个向量
    叉积判断旋转方向，确保返回正值角度
    """
    norm = (qx * qx + qy * qy) * (rx * rx + ry * ry)
    if norm <= 0.0:
        return 0.0
    cos_angle = (qx * rx + qy * ry) / math.sqrt(norm)
    if cos_angle >= 1.0:
        return 0.0
    if cos_angle <= -1.0:
        return math.pi
    value = math.acos(cos_angle)
    if qx * ry < qy * rx:
        value = -value
    if value < 0.0:
        value += TWO_PI
    return value


def angle_points(p: Point2D, q: Point2D, r: Point2D) -> float:
    """计算点 q 相对于点 p 的方向到点 r 相对于点 p 的方向的夹角"""
    return angle_vectors(q.x - p.x, q.y - p.y, r.x - p.x, r.y - p.y)


def same_point(a: Point2D, b: Point2D) -> bool:
    """判断两个点是否（近似）相同"""
    return abs(a.x - b.x) < EPSILON and abs(a.y - b.y) < EPSILON


def log2(value: float) -> float:
    """计算以 2 为底的对数"""
    return math.log(value) / math.log(2.0)


def f_attr_scalar(distance: float, ideal_edge_length: float) -> float:
    """
    计算吸引力标量值（力导向布局中的弹簧力）
    - distance: 当前两点间的实际距离
    - ideal_edge_length: 理想边长度
    基于对数函数的吸引力模型，距离越远吸引力越大
    """
    if distance <= 0.0:
        return -1e10
    c = log2(distance / ideal_edge_length)
    return c * distance * distance / (ideal_edge_length * ideal_edge_length * ideal_edge_length)


def get_post_rep_force_strength(node_count: int) -> float:
    """
    获取后处理阶段的斥力强度
    节点数越多，斥力强度越小（避免过度推散），上限 0.2
    """
    return min(0.2, 400.0 / max(node_count, 1))


def max_radius(box_length: float, iteration: int) -> float:
    """
    限制节点单次移动的最大半径
    - 第 1 次迭代：使用很小的半径（box_length/1000），避免初始布局剧烈震荡
    - 后续迭代：使用较大半径（box_length/5），允许更大的调整幅度
    """
    return box_length / 1000.0 if iteration == 1 else box_length / 5.0


def edge_key(a_id: int, b_id: int) -> tuple[int, int]:
    """生成无向边的键（较小的 id 在前），用于边的去重和查找"""
    return (min(a_id, b_id), max(a_id, b_id))


def segment_intersects(a1: Point2D, a2: Point2D, b1: Point2D, b2: Point2D) -> bool:
    """
    判断两条线段是否相交
    使用向量叉积方法判断方向，用于 fix_twisted_splits 中检测交叉边
    """
    def orientation(p: Point2D, q: Point2D, r: Point2D) -> float:
        """计算三点 p→q→r 的方向（正=逆时针，负=顺时针，零=共线）"""
        return (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)

    def on_segment(p: Point2D, q: Point2D, r: Point2D) -> bool:
        """判断点 q 是否在线段 pr 上"""
        return min(p.x, r.x) - EPSILON <= q.x <= max(p.x, r.x) + EPSILON and min(p.y, r.y) - EPSILON <= q.y <= max(p.y, r.y) + EPSILON

    o1 = orientation(a1, a2, b1)
    o2 = orientation(a1, a2, b2)
    o3 = orientation(b1, b2, a1)
    o4 = orientation(b1, b2, a2)
    if (o1 > 0 and o2 < 0 or o1 < 0 and o2 > 0) and (o3 > 0 and o4 < 0 or o3 < 0 and o4 > 0):
        return True
    if abs(o1) < EPSILON and on_segment(a1, b1, a2):
        return True
    if abs(o2) < EPSILON and on_segment(a1, b2, a2):
        return True
    if abs(o3) < EPSILON and on_segment(b1, a1, b2):
        return True
    if abs(o4) < EPSILON and on_segment(b1, a2, b2):
        return True
    return False


# ==================== 渲染选项枚举 ====================

class NodeLabelMode(Enum):
    """节点标签显示模式"""
    NONE = "none"            # 不显示标签
    NAME = "name"            # 仅显示节点名称
    DEPTH = "depth"          # 仅显示深度值
    LENGTH = "length"        # 仅显示节点长度
    NAME_DEPTH = "name-depth"  # 显示名称和深度（默认）
    NAME_LENGTH = "name-length"  # 显示名称和长度
    NAME_DEPTH_LENGTH = "name-depth-length"  # 显示名称、深度和长度


class EdgeLabelMode(Enum):
    """边标签显示模式"""
    NONE = "none"          # 不显示标签
    NAME = "name"          # 显示边名称（默认）
    DERIVED = "derived"    # 显示派生的连接描述（如 "1+ -> 2+"）


@dataclass(frozen=True)
class BandageRenderOptions:
    """渲染选项配置，控制节点/边标签的显示方式和字体"""
    node_label_mode: NodeLabelMode = NodeLabelMode.NAME_DEPTH  # 节点标签模式
    edge_label_mode: EdgeLabelMode = EdgeLabelMode.NAME        # 边标签模式
    node_font_size: int = 14   # 节点标签字号
    edge_font_size: int = 12   # 边标签字号
    font_color: RgbaColor = BLACK  # 标签颜色

    @staticmethod
    def parse_node_label_mode(value: str) -> NodeLabelMode:
        """从字符串解析节点标签模式（支持多种别名）"""
        normalized = value.lower()
        if normalized in {"none", "off", "hide"}:
            return NodeLabelMode.NONE
        if normalized in {"name", "id"}:
            return NodeLabelMode.NAME
        if normalized in {"depth", "dp"}:
            return NodeLabelMode.DEPTH
        if normalized in {"length", "len"}:
            return NodeLabelMode.LENGTH
        if normalized in {"name-depth", "depth-name", "all", "full"}:
            return NodeLabelMode.NAME_DEPTH
        if normalized in {"name-length", "length-name"}:
            return NodeLabelMode.NAME_LENGTH
        if normalized in {"name-depth-length", "name-length-depth", "depth-length-name", "full-detail"}:
            return NodeLabelMode.NAME_DEPTH_LENGTH
        raise ValueError(f"Unsupported node label mode: {value}. Expected none|name|depth|length|name-depth|name-length|name-depth-length")

    @staticmethod
    def parse_edge_label_mode(value: str) -> EdgeLabelMode:
        """从字符串解析边标签模式（支持多种别名）"""
        normalized = value.lower()
        if normalized in {"none", "off", "hide"}:
            return EdgeLabelMode.NONE
        if normalized in {"name", "id", "on", "show"}:
            return EdgeLabelMode.NAME
        if normalized in {"derived", "connection", "connections"}:
            return EdgeLabelMode.DERIVED
        raise ValueError(f"Unsupported edge label mode: {value}. Expected none|name|derived")


# ==================== 图布局设置 ====================

@dataclass
class BandageSettings:
    """
    Bandage 渲染和布局的所有可配置参数
    包括布局算法参数、节点/边的视觉样式参数等
    """
    double_mode: bool = False               # 双链模式：同时显示正链(+)和负链(-)节点
    arrowheads_in_single_mode: bool = False # 单链模式下是否显示箭头
    mean_node_length: float = 40.0          # 节点平均绘制长度（像素）
    min_total_graph_length: float = 500.0   # 整图最小绘制长度
    minimum_node_length: float = 5.0        # 单个节点最小绘制长度
    edge_length: float = 5.0                # 边的默认长度
    node_segment_length: float = 20.0       # 节点被分割为多段布局时，每段的长度
    component_separation: float = 50.0      # 连通分量之间的间距
    average_node_width: float = 5.0         # 节点平均宽度（与深度相关）
    depth_effect_on_width: float = 0.5      # 深度对节点宽度的影响系数
    depth_power: float = 0.5                # 深度幂次，控制宽度随深度的非线性变化
    edge_width: float = 1.5                 # 边的线宽
    outline_thickness: float = 0.0          # 节点轮廓粗细
    selection_thickness: float = 1.0        # 选中元素的轮廓粗细
    double_mode_node_separation: float = 2.0  # 双链模式下正负链节点的间距
    graph_layout_quality: int = 2           # 布局质量等级（0=低, 4=高，默认2=中等）
    linear_layout: bool = False             # 是否使用线性布局（节点排成一行）
    use_auto_node_length: bool = True       # 是否自动计算节点长度（基于图的总长度）
    manual_node_length_per_megabase: float = 10000.0  # 手动模式下每 Mb 对应的像素长度
    auto_node_length_per_megabase: float = 10000.0    # 自动模式下每 Mb 对应的像素长度（会根据图大小自动调整）
    random_colour_positive_saturation: int = 127  # 随机颜色的饱和度
    random_colour_positive_lightness: int = 150    # 随机颜色的亮度
    random_colour_positive_opacity: int = 255      # 随机颜色的不透明度
    edge_colour: RgbaColor = field(default_factory=lambda: RgbaColor(0, 0, 0, 180))      # 边颜色（半透明黑色）
    outline_colour: RgbaColor = field(default_factory=lambda: RgbaColor(0, 0, 0, 255))   # 轮廓颜色
    selection_colour: RgbaColor = field(default_factory=lambda: RgbaColor(0, 0, 255, 255)) # 选中高亮颜色
    background_colour: RgbaColor = field(default_factory=lambda: WHITE)                    # 背景颜色

    def has_arrow_in_single_mode(self) -> bool:
        """判断单链模式下是否需要显示箭头"""
        return self.double_mode or self.arrowheads_in_single_mode


# ==================== 布局图数据结构 ====================

class LayoutNode:
    """
    布局图中的节点
    用于力导向布局算法，存储节点的位置、受力等信息
    - x, y: 节点当前坐标
    - force_x, force_y: 当前迭代中的合力分量
    - last_move_x, last_move_y: 上一次迭代的位移（用于防振荡）
    - attr_x, attr_y: 吸引力分量
    - rep_x, rep_y: 斥力分量
    """
    def __init__(self, node_id: int) -> None:
        self.id = node_id
        self.edges: list[LayoutEdge] = []
        self.x = 0.0
        self.y = 0.0
        self.force_x = 0.0
        self.force_y = 0.0
        self.last_move_x = 0.0
        self.last_move_y = 0.0
        self.attr_x = 0.0
        self.attr_y = 0.0
        self.rep_x = 0.0
        self.rep_y = 0.0


class LayoutEdge:
    """
    布局图中的边
    - source, target: 边的起止节点
    - desired_length: 期望的边长度（布局目标）
    - internal: 是否为节点内部边（将长节点分割为多段时的连接边）
    """
    def __init__(self, source: LayoutNode, target: LayoutNode, desired_length: float, internal: bool) -> None:
        self.source = source
        self.target = target
        self.desired_length = desired_length
        self.internal = internal


class LayoutGraph:
    """
    布局图，由 LayoutNode 和 LayoutEdge 组成
    是力导向布局算法的输入图结构
    """
    def __init__(self) -> None:
        self.nodes: list[LayoutNode] = []
        self.edges: list[LayoutEdge] = []

    def new_node(self) -> LayoutNode:
        """创建新节点并添加到图中"""
        node = LayoutNode(len(self.nodes))
        self.nodes.append(node)
        return node

    def new_edge(self, a: LayoutNode, b: LayoutNode, desired_length: float, internal: bool) -> LayoutEdge:
        """创建新边并添加到图中，同时在起止节点上记录该边"""
        edge = LayoutEdge(a, b, desired_length, internal)
        self.edges.append(edge)
        a.edges.append(edge)
        b.edges.append(edge)
        return edge

    def clear(self) -> None:
        """清空图中所有节点和边"""
        self.nodes.clear()
        self.edges.clear()


# ==================== DeBruijn 图数据结构 ====================

class OgdfNode:
    """
    OGDF 节点包装器
    一个 DeBruijn 图中的节点在布局图中可能被分割为多个 LayoutNode
    （因为长节点需要多段来表示其弯曲路径），OgdfNode 管理这些分割后的节点序列
    """
    def __init__(self) -> None:
        self._nodes: list[LayoutNode] = []

    def add_ogdf_node(self, node: LayoutNode) -> None:
        """添加一个布局节点到序列中"""
        self._nodes.append(node)

    def get_first(self) -> LayoutNode | None:
        """获取序列中第一个节点（节点路径的起点）"""
        return self._nodes[0] if self._nodes else None

    def get_second(self) -> LayoutNode | None:
        """获取序列中第二个节点；若只有1个则返回第一个"""
        return self._nodes[1] if len(self._nodes) >= 2 else self.get_first()

    def get_last(self) -> LayoutNode | None:
        """获取序列中最后一个节点（节点路径的终点）"""
        return self._nodes[-1] if self._nodes else None

    def get_second_last(self) -> LayoutNode | None:
        """获取序列中倒数第二个节点；若只有1个则返回最后一个"""
        return self._nodes[-2] if len(self._nodes) >= 2 else self.get_last()

    def nodes(self) -> list[LayoutNode]:
        """返回所有布局节点"""
        return self._nodes


class DeBruijnNode:
    """
    DeBruijn 图中的节点
    代表 GFA 文件中的一个序列片段（S 行），带有方向标识（+/-）
    - name: 节点名称（如 "1+", "2-"）
    - depth: 测序深度（覆盖度）
    - sequence: DNA 序列
    - length: 序列长度
    - reverse_complement: 反向互补节点的引用（如 1+ 的反向互补是 1-）
    - ogdf_node: 对应的布局节点序列
    - drawn: 是否被标记为需要绘制
    """
    def __init__(self, name: str, depth: float, sequence: str, length: int) -> None:
        self.name = name
        self.depth = depth
        self.sequence = sequence or ""
        self.length = max(length, 0)
        self.edges: list[DeBruijnEdge] = []
        self.depth_relative_to_mean_drawn_depth = 1.0  # 相对于平均绘制深度的比率
        self.width = 0.0        # 节点绘制宽度（基于深度计算）
        self.color = RgbaColor(128, 128, 128, 255)  # 节点颜色
        self.reverse_complement: DeBruijnNode | None = None  # 反向互补节点
        self.ogdf_node: OgdfNode | None = None  # 对应的布局节点序列
        self.drawn = False       # 是否标记为需要绘制

    def get_name_without_sign(self) -> str:
        """获取不带方向标记的节点名称（如 "1+" → "1"）"""
        return self.name[:-1] if self.name else ""

    def add_edge(self, edge: DeBruijnEdge) -> None:
        """将与该节点相连的边添加到边列表（去重）"""
        if edge not in self.edges:
            self.edges.append(edge)

    def get_upstream_nodes(self) -> list[DeBruijnNode]:
        """获取所有上游节点（即所有指向该节点的边的起始节点）"""
        upstream_nodes: list[DeBruijnNode] = []
        for edge in self.edges:
            if self is edge.ending_node:
                upstream_nodes.append(edge.starting_node)
        return upstream_nodes

    def is_positive_node(self) -> bool:
        """判断是否为正链节点（名称以 '+' 结尾）"""
        return bool(self.name) and self.name[-1] == "+"

    def is_negative_node(self) -> bool:
        """判断是否为负链节点（名称以 '-' 结尾）"""
        return bool(self.name) and self.name[-1] == "-"

    def in_ogdf(self) -> bool:
        """判断该节点是否已添加到布局图中"""
        return self.ogdf_node is not None

    def this_or_reverse_complement_in_ogdf(self) -> bool:
        """判断该节点或其反向互补节点是否已在布局图中"""
        return self.in_ogdf() or (self.reverse_complement is not None and self.reverse_complement.in_ogdf())

    def this_or_reverse_complement_not_in_ogdf(self) -> bool:
        """判断该节点和其反向互补节点是否都不在布局图中"""
        return not self.this_or_reverse_complement_in_ogdf()

    def set_as_drawn(self) -> None:
        """标记该节点为需要绘制"""
        self.drawn = True

    def reset_node(self) -> None:
        """重置节点的布局状态（用于重新布局）"""
        self.ogdf_node = None
        self.drawn = False

    def add_to_ogdf_graph(self, ogdf_graph: LayoutGraph, x_pos: float, y_pos: float, settings: BandageSettings) -> None:
        """
        将该 DeBruijn 节点转换为布局图中的节点和边
        长节点会被分割为多个段（每段 node_segment_length），段之间用内部边连接
        - x_pos, y_pos: 线性布局时的起始坐标
        """
        if self.this_or_reverse_complement_in_ogdf():
            return  # 避免重复添加（正链和负链共享布局位置）
        self.ogdf_node = OgdfNode()
        drawn_node_length = self.get_drawn_node_length(settings)
        number_of_graph_edges = self.get_number_of_ogdf_graph_edges(drawn_node_length, settings)
        number_of_graph_nodes = number_of_graph_edges + 1
        drawn_length_per_edge = drawn_node_length / number_of_graph_edges
        previous_node: LayoutNode | None = None
        for _ in range(number_of_graph_nodes):
            new_node = ogdf_graph.new_node()
            if settings.linear_layout:
                # 线性布局：节点沿 x 轴依次排列
                new_node.x = x_pos
                new_node.y = y_pos
                x_pos += settings.node_segment_length
            self.ogdf_node.add_ogdf_node(new_node)
            if previous_node is not None:
                # 在相邻段之间添加内部边
                ogdf_graph.new_edge(previous_node, new_node, drawn_length_per_edge, True)
            previous_node = new_node

    def get_drawn_node_length(self, settings: BandageSettings) -> float:
        """计算该节点在图中的绘制长度（基于序列长度和比例因子）"""
        drawn_node_length = self.get_node_length_per_megabase(settings) * self.length / 1_000_000.0
        return max(drawn_node_length, settings.minimum_node_length)

    def get_number_of_ogdf_graph_edges(self, drawn_node_length: float, settings: BandageSettings) -> int:
        """根据绘制长度和段长度，计算需要分割成多少条内部边"""
        return max(int(math.ceil(drawn_node_length / settings.node_segment_length)), 1)

    def get_node_length_per_megabase(self, settings: BandageSettings) -> float:
        """获取每 Mb 序列对应的像素长度（自动或手动模式）"""
        return settings.auto_node_length_per_megabase if settings.use_auto_node_length else settings.manual_node_length_per_megabase

    def sequence_is_missing(self) -> bool:
        """判断序列数据是否缺失"""
        return not self.sequence

    def build_display_points(self) -> list[Point2D]:
        """
        构建该节点在画布上的显示路径点序列
        如果该节点在布局图中，直接使用布局节点坐标
        如果只有反向互补节点在布局图中，则反转使用其坐标
        """
        line_points: list[Point2D] = []
        if self.ogdf_node is not None:
            for node in self.ogdf_node.nodes():
                line_points.append(Point2D(node.x, node.y))
            return line_points
        if self.reverse_complement is None or self.reverse_complement.ogdf_node is None:
            return line_points
        # 反向互补节点在布局图中，反转坐标序列
        for node in reversed(self.reverse_complement.ogdf_node.nodes()):
            line_points.append(Point2D(node.x, node.y))
        return line_points

    @staticmethod
    def get_node_width(depth_relative_to_mean_drawn_depth: float, settings: BandageSettings) -> float:
        """
        根据深度计算节点绘制宽度
        深度越大，节点越宽（使用幂函数实现非线性关系）
        """
        clamped_depth = max(0.0, depth_relative_to_mean_drawn_depth)
        width_relative_to_average = (math.pow(clamped_depth, settings.depth_power) - 1.0) * settings.depth_effect_on_width + 1.0
        return settings.average_node_width * width_relative_to_average


class DeBruijnEdge:
    """
    DeBruijn 图中的边
    代表 GFA 文件中的一个连接（L 行），连接两个有方向的节点
    - starting_node, ending_node: 边的起止节点
    - overlap: 重叠长度（CIGAR 解析得到）
    - label: 边的标签（来自 GFA 的 ID/Name/LB 等标签）
    - reverse_complement: 反向互补边的引用
    - drawn: 是否需要绘制
    """
    def __init__(self, starting_node: DeBruijnNode, ending_node: DeBruijnNode, overlap: int, label: str | None) -> None:
        self.starting_node = starting_node
        self.ending_node = ending_node
        self.overlap = overlap
        self.label = label
        self.reverse_complement: DeBruijnEdge | None = None
        self.drawn = False

    def reset(self) -> None:
        """重置边的绘制状态"""
        self.drawn = False

    def determine_if_drawn(self, settings: BandageSettings) -> None:
        """根据设置判断该边是否需要绘制"""
        self.drawn = self.edge_is_visible(settings)

    def is_positive_edge(self) -> bool:
        """
        判断是否为"正向边"
        在单链模式下，只绘制正向边以避免重复
        规则：正链→正链 的边为正向；负链→负链 的边为反向
        混合方向边通过节点名称大小比较决定
        """
        if self.starting_node.is_positive_node() and self.ending_node.is_positive_node():
            return True
        if self.starting_node.is_negative_node() and self.ending_node.is_negative_node():
            return False
        if self.is_own_reverse_complement():
            return True
        if self.reverse_complement is None:
            return True
        return self.starting_node.name > self.reverse_complement.starting_node.name

    def is_own_reverse_complement(self) -> bool:
        """判断该边是否是其自身的反向互补（自环边的情况）"""
        return self is self.reverse_complement

    def add_to_ogdf_graph(self, ogdf_graph: LayoutGraph, settings: BandageSettings) -> None:
        """
        将该边添加到布局图中
        确定起止节点在布局图中的对应端点，然后添加外部边
        正链节点的边从末尾出发，负链节点的边从头部出发
        """
        # 确定起始端在布局图中的位置
        if self.starting_node.in_ogdf():
            first_edge_ogdf_node = self.starting_node.ogdf_node.get_last() if self.starting_node.ogdf_node else None
        elif self.starting_node.reverse_complement is not None and self.starting_node.reverse_complement.in_ogdf():
            first_edge_ogdf_node = self.starting_node.reverse_complement.ogdf_node.get_first() if self.starting_node.reverse_complement.ogdf_node else None
        else:
            return
        # 确定终止端在布局图中的位置
        if self.ending_node.in_ogdf():
            second_edge_ogdf_node = self.ending_node.ogdf_node.get_first() if self.ending_node.ogdf_node else None
        elif self.ending_node.reverse_complement is not None and self.ending_node.reverse_complement.in_ogdf():
            second_edge_ogdf_node = self.ending_node.reverse_complement.ogdf_node.get_last() if self.ending_node.reverse_complement.ogdf_node else None
        else:
            return
        if first_edge_ogdf_node is None or second_edge_ogdf_node is None:
            return
        # 处理自环边：如果节点只有1条内部边，无法画出回环
        if self.starting_node is self.ending_node:
            node_graph_edges = self.starting_node.get_number_of_ogdf_graph_edges(self.starting_node.get_drawn_node_length(settings), settings)
            if node_graph_edges == 1:
                return
        ogdf_graph.new_edge(first_edge_ogdf_node, second_edge_ogdf_node, settings.edge_length, False)

    def edge_is_visible(self, settings: BandageSettings) -> bool:
        """
        判断该边是否应该可见
        - 双链模式：两端节点都标记绘制时可见
        - 单链模式：还需满足是正向边
        """
        if settings.double_mode:
            return self.starting_node.drawn and self.ending_node.drawn
        draw_edge = (
            (self.starting_node.drawn or (self.starting_node.reverse_complement is not None and self.starting_node.reverse_complement.drawn))
            and (self.ending_node.drawn or (self.ending_node.reverse_complement is not None and self.ending_node.reverse_complement.drawn))
        )
        return draw_edge and self.is_positive_edge()


# ==================== 组装图（核心图结构） ====================

class AssemblyGraph:
    """
    组装图：整个 DeBruijn 图的管理类
    负责从 GFA 文件解析图结构、构建布局图、执行布局计算等
    """
    def __init__(self, settings: BandageSettings, rng: random.Random | None = None) -> None:
        self.settings = settings
        self.de_bruijn_graph_nodes: dict[str, DeBruijnNode] = {}   # 所有节点（按名称索引）
        self.de_bruijn_graph_edges: dict[tuple[str, str], DeBruijnEdge] = {}  # 所有边（按起止节点名索引）
        self.ogdf_graph = LayoutGraph()  # 用于布局计算的图
        self.node_count = 0    # 正链节点数
        self.edge_count = 0    # 正向边数
        self.total_length = 0  # 正链节点总长度
        self.mean_depth = 0.0  # 加权平均深度
        self.random = rng if rng is not None else random.Random(time.time_ns())

    @classmethod
    def load_from_gfa(cls, input_path: Path, settings: BandageSettings, rng: random.Random | None = None) -> "AssemblyGraph":
        """从 GFA 文件加载组装图"""
        graph = cls(settings, rng)
        graph.build_de_bruijn_graph_from_gfa(input_path)
        graph.determine_graph_info()
        return graph

    def prepare_for_rendering(self) -> None:
        """准备渲染：计算图信息、深度相对值、节点宽度和颜色"""
        self.determine_graph_info()
        self.recalculate_all_depths_relative_to_drawn_mean()
        self.recalculate_all_node_widths()
        self.reset_all_node_colours()

    def render_whole_graph(self) -> None:
        """渲染整张图：清空旧布局 → 构建布局图 → 计算布局 → 更新样式"""
        self.clear_ogdf_graph_and_reset_nodes()
        self.build_ogdf_graph_from_nodes_and_edges()
        self.layout_graph()
        self.recalculate_all_depths_relative_to_drawn_mean()
        self.recalculate_all_node_widths()
        self.reset_all_node_colours()

    def get_displayed_nodes(self) -> list[DeBruijnNode]:
        """获取所有标记为需要绘制的节点"""
        return [node for node in self.de_bruijn_graph_nodes.values() if node.drawn]

    def get_displayed_edges(self) -> list[DeBruijnEdge]:
        """获取所有标记为需要绘制的边"""
        return [edge for edge in self.de_bruijn_graph_edges.values() if edge.drawn]

    def build_de_bruijn_graph_from_gfa(self, input_path: Path) -> None:
        """
        从 GFA 文件构建 DeBruijn 图
        解析 GFA 格式：
        - S 行（Segment）：节点，包含名称、序列、长度(LN)、深度(DP/KC/RC/FC)等标签
        - L 行（Link）：边，包含起始节点+方向、终止节点+方向、CIGAR 重叠信息
        同时自动创建缺失的反向互补节点
        """
        edge_starting_node_names: list[str] = []
        edge_ending_node_names: list[str] = []
        edge_overlaps: list[int] = []
        edge_labels: list[str | None] = []

        with input_path.open("r", encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\n")
                if not line:
                    continue
                line_parts = line.split("\t")
                if not line_parts:
                    continue
                record_type = line_parts[0]
                if record_type == "S":
                    # 解析 Segment 行：S  节点名  序列  [可选标签...]
                    if len(line_parts) < 3:
                        continue
                    node_name = line_parts[1]
                    last_char = node_name[-1] if node_name else "+"
                    if last_char not in {"+", "-"}:
                        node_name += "+"  # GFA1 格式默认为正链
                    sequence = line_parts[2]
                    # 解析各种可选标签
                    ln = parse_int_tag(line_parts, "LN:i:", parse_int_tag(line_parts, "ln:i:", 0))
                    dp = parse_float_tag(
                        line_parts,
                        "DP:f:",
                        parse_float_tag(line_parts, "dp:f:", parse_float_tag(line_parts, "dp:i:", math.nan)),
                    )
                    kc = parse_int_tag(line_parts, "KC:i:", parse_int_tag(line_parts, "kc:i:", -1))
                    rc = parse_int_tag(line_parts, "RC:i:", parse_int_tag(line_parts, "rc:i:", -1))
                    fc = parse_int_tag(line_parts, "FC:i:", parse_int_tag(line_parts, "fc:i:", -1))
                    if sequence == "*" or not sequence:
                        # 序列缺失时使用 LN 标签的长度
                        length = ln
                        sequence = ""
                    else:
                        length = len(sequence)
                    # 计算节点深度：优先使用 DP 标签，其次从 KC/RC/FC 标签推算
                    node_depth = 1.0
                    if not math.isnan(dp):
                        node_depth = dp
                    elif kc >= 0 and length > 0:
                        node_depth = kc / length
                    elif rc >= 0 and length > 0:
                        node_depth = rc / length
                    elif fc >= 0 and length > 0:
                        node_depth = fc / length
                    self.de_bruijn_graph_nodes[node_name] = DeBruijnNode(node_name, node_depth, sequence, length)
                elif record_type == "L":
                    # 解析 Link 行：L  起始名  起始方向  终止名  终止方向  CIGAR  [可选标签...]
                    if len(line_parts) < 6:
                        continue
                    starting_node = line_parts[1] + line_parts[2]  # 如 "1+"
                    ending_node = line_parts[3] + line_parts[4]    # 如 "2-"
                    edge_starting_node_names.append(starting_node)
                    edge_ending_node_names.append(ending_node)
                    edge_labels.append(parse_string_tag(line_parts, "ID:Z:", "id:Z:", "Name:Z:", "NAME:Z:", "name:Z:", "LB:Z:", "lb:Z:"))
                    cigar = line_parts[5]
                    # 解析 CIGAR 字符串获取重叠长度
                    if cigar == "*":
                        edge_overlaps.append(0)
                    elif cigar_contains_only_m(cigar):
                        edge_overlaps.append(get_length_from_simple_cigar(cigar))
                    else:
                        edge_overlaps.append(get_length_from_cigar(cigar))

        # 为每个节点创建反向互补节点（如果 GFA 中未显式定义）
        for node in list(self.de_bruijn_graph_nodes.values()):
            self.make_reverse_complement_node_if_necessary(node)
        # 建立正负链节点之间的互相引用
        self.point_each_node_to_its_reverse_complement()

        # 创建所有边（包括正向和反向互补边）
        for index, starting_name in enumerate(edge_starting_node_names):
            self.create_de_bruijn_edge(starting_name, edge_ending_node_names[index], edge_overlaps[index], edge_labels[index])

    def create_de_bruijn_edge(self, node1_name: str, node2_name: str, overlap: int, label: str | None) -> None:
        """
        创建一条 DeBruijn 边及其反向互补边
        正向边连接 node1 → node2，反向互补边连接 node2的反向 → node1的反向
        如果正向边和反向互补边是同一条边（自反向互补），则只创建一条
        """
        node1_opposite = get_opposite_node_name(node1_name)
        node2_opposite = get_opposite_node_name(node2_name)
        node1 = self.de_bruijn_graph_nodes.get(node1_name)
        node2 = self.de_bruijn_graph_nodes.get(node2_name)
        neg_node1 = self.de_bruijn_graph_nodes.get(node1_opposite)
        neg_node2 = self.de_bruijn_graph_nodes.get(node2_opposite)
        if node1 is None or node2 is None or neg_node1 is None or neg_node2 is None:
            return
        existing_key = (node1.name, node2.name)
        if existing_key in self.de_bruijn_graph_edges:
            return  # 避免重复边
        is_own_pair = node1 is neg_node2 and node2 is neg_node1  # 正反向是同一个边
        forward_edge = DeBruijnEdge(node1, node2, overlap, label)
        backward_edge = forward_edge if is_own_pair else DeBruijnEdge(neg_node2, neg_node1, overlap, label)
        forward_edge.reverse_complement = backward_edge
        backward_edge.reverse_complement = forward_edge
        self.de_bruijn_graph_edges[(forward_edge.starting_node.name, forward_edge.ending_node.name)] = forward_edge
        if not is_own_pair:
            self.de_bruijn_graph_edges[(backward_edge.starting_node.name, backward_edge.ending_node.name)] = backward_edge
        node1.add_edge(forward_edge)
        node2.add_edge(forward_edge)
        neg_node1.add_edge(backward_edge)
        neg_node2.add_edge(backward_edge)

    def clear_ogdf_graph_and_reset_nodes(self) -> None:
        """清空布局图并重置所有节点的布局状态"""
        for node in self.de_bruijn_graph_nodes.values():
            node.reset_node()
        for edge in self.de_bruijn_graph_edges.values():
            edge.reset()
        self.ogdf_graph.clear()

    def determine_graph_info(self) -> None:
        """
        计算图的基本统计信息
        - 正链节点数、正向边数、正链总长度
        - 长度加权平均深度
        - 自动计算节点长度比例（使整图适配目标绘制长度）
        """
        positive_node_count = 0
        positive_total_length = 0
        weighted_length = 0
        weighted_depth = 0.0
        for node in self.de_bruijn_graph_nodes.values():
            node_length = node.length
            weighted_length += node_length
            weighted_depth += node_length * node.depth
            if node.is_positive_node():
                positive_total_length += node_length
                positive_node_count += 1
        positive_edge_count = 0
        for edge in self.de_bruijn_graph_edges.values():
            if edge.is_positive_edge():
                positive_edge_count += 1
        self.node_count = positive_node_count
        self.edge_count = positive_edge_count
        self.total_length = positive_total_length
        self.mean_depth = 0.0 if weighted_length == 0 else weighted_depth / weighted_length
        # 根据目标绘制长度和总序列长度，计算每 Mb 对应的像素数
        target_drawn_graph_length = max(self.node_count * self.settings.mean_node_length, self.settings.min_total_graph_length)
        megabases = positive_total_length / 1_000_000.0
        self.settings.auto_node_length_per_megabase = target_drawn_graph_length / megabases if megabases > 0.0 else 10_000.0

    def build_ogdf_graph_from_nodes_and_edges(self) -> None:
        """
        将 DeBruijn 图转换为布局图
        - 在单链模式下只标记正链节点为 drawn
        - 将标记为 drawn 的节点添加到布局图中
        - 将标记为 drawn 的边添加到布局图中
        """
        for node in self.de_bruijn_graph_nodes.values():
            if node.is_positive_node() or self.settings.double_mode:
                node.set_as_drawn()

        if self.use_linear_layout():
            # 线性布局：按节点名排序，沿 x 轴排列
            sorted_drawn_nodes: list[DeBruijnNode] = []
            numeric_sort = True
            for node in self.de_bruijn_graph_nodes.values():
                if node.drawn and node.this_or_reverse_complement_not_in_ogdf():
                    try:
                        int(node.get_name_without_sign())
                    except ValueError:
                        numeric_sort = False
                        break
                    sorted_drawn_nodes.append(node)
            if numeric_sort:
                sorted_drawn_nodes.sort(key=lambda node: int(node.get_name_without_sign()))
            else:
                sorted_drawn_nodes = [
                    node
                    for node in self.de_bruijn_graph_nodes.values()
                    if node.drawn and node.this_or_reverse_complement_not_in_ogdf()
                ]
                sorted_drawn_nodes.sort(key=lambda node: node.get_name_without_sign().upper())

            last_x_pos = 0.0
            for node in sorted_drawn_nodes:
                if node.this_or_reverse_complement_in_ogdf():
                    continue
                # 根据上游节点的位置确定当前节点的起始 x 坐标
                for upstream_node in node.get_upstream_nodes():
                    if not upstream_node.in_ogdf() or upstream_node.ogdf_node is None or upstream_node.ogdf_node.get_last() is None:
                        continue
                    last_x_pos = max(last_x_pos, upstream_node.ogdf_node.get_last().x)
                x_pos = last_x_pos + self.settings.edge_length
                node.add_to_ogdf_graph(self.ogdf_graph, x_pos, 0.0, self.settings)
                if node.ogdf_node is not None and node.ogdf_node.get_last() is not None:
                    last_x_pos = node.ogdf_node.get_last().x
        else:
            # 力导向布局：所有节点初始位置为 (0, 0)
            for node in self.de_bruijn_graph_nodes.values():
                if node.drawn and node.this_or_reverse_complement_not_in_ogdf():
                    node.add_to_ogdf_graph(self.ogdf_graph, 0.0, 0.0, self.settings)

        for edge in self.de_bruijn_graph_edges.values():
            edge.determine_if_drawn(self.settings)
            if edge.drawn:
                edge.add_to_ogdf_graph(self.ogdf_graph, self.settings)

    def layout_graph(self) -> None:
        """使用力导向布局算法计算节点位置"""
        worker = GraphLayoutWorker(
            self.ogdf_graph,
            self.settings.graph_layout_quality,
            self.use_linear_layout(),
            self.settings.component_separation,
            rng=self.random,
        )
        worker.layout_graph()

    def recalculate_all_depths_relative_to_drawn_mean(self) -> None:
        """重新计算所有节点相对于平均绘制深度的比率（用于确定节点宽度）"""
        mean_drawn_depth = self.get_mean_depth(True)
        for node in self.de_bruijn_graph_nodes.values():
            node.depth_relative_to_mean_drawn_depth = 1.0 if mean_drawn_depth == 0.0 else node.depth / mean_drawn_depth

    def recalculate_all_node_widths(self) -> None:
        """根据深度重新计算所有节点的绘制宽度"""
        for node in self.de_bruijn_graph_nodes.values():
            node.width = DeBruijnNode.get_node_width(node.depth_relative_to_mean_drawn_depth, self.settings)

    def reset_all_node_colours(self) -> None:
        """为所有绘制节点分配随机颜色（正链和负链节点共享同一颜色）"""
        for node in self.de_bruijn_graph_nodes.values():
            if not node.drawn:
                continue
            hue = self.random.randrange(360)
            positive_colour = hsl_to_rgb(
                hue / 360.0,
                self.settings.random_colour_positive_saturation / 255.0,
                self.settings.random_colour_positive_lightness / 255.0,
                self.settings.random_colour_positive_opacity,
            )
            node.color = positive_colour
            if node.reverse_complement is not None:
                node.reverse_complement.color = positive_colour

    def use_linear_layout(self) -> bool:
        """判断是否使用线性布局（无边或用户指定线性布局时）"""
        return self.edge_count == 0 or self.settings.linear_layout

    def get_mean_depth(self, drawn_nodes_only: bool) -> float:
        """计算长度加权平均深度"""
        total_node_length = 0
        depth_sum = 0.0
        for node in self.de_bruijn_graph_nodes.values():
            if drawn_nodes_only and not node.drawn:
                continue
            total_node_length += node.length
            depth_sum += node.length * node.depth
        return 0.0 if total_node_length == 0 else depth_sum / total_node_length

    def make_reverse_complement_node_if_necessary(self, node: DeBruijnNode) -> None:
        """如果反向互补节点不存在，则自动创建（计算反向互补序列）"""
        reverse_complement_name = get_opposite_node_name(node.name)
        if reverse_complement_name in self.de_bruijn_graph_nodes:
            return
        reverse_complement_sequence = "" if node.sequence_is_missing() else reverse_complement(node.sequence)
        self.de_bruijn_graph_nodes[reverse_complement_name] = DeBruijnNode(reverse_complement_name, node.depth, reverse_complement_sequence, node.length)

    def point_each_node_to_its_reverse_complement(self) -> None:
        """建立每个正链节点与其负链节点之间的双向引用"""
        for node in self.de_bruijn_graph_nodes.values():
            if not node.is_positive_node():
                continue
            reverse_node = self.de_bruijn_graph_nodes.get(get_opposite_node_name(node.name))
            if reverse_node is not None:
                node.reverse_complement = reverse_node
                reverse_node.reverse_complement = node


# ==================== GFA 解析辅助函数 ====================

def get_opposite_node_name(node_name: str | None) -> str:
    """获取节点的反向互补节点名称（切换 + / - 方向）"""
    if not node_name:
        return "+"
    last_char = node_name[-1]
    if last_char == "+":
        return node_name[:-1] + "-"
    if last_char == "-":
        return node_name[:-1] + "+"
    return node_name + "-"


def reverse_complement(forward_sequence: str) -> str:
    """计算 DNA 序列的反向互补（A↔T, G↔C, 其他→N）"""
    output: list[str] = []
    for letter in reversed(forward_sequence):
        upper = letter.upper()
        if upper == "A":
            output.append("T")
        elif upper == "T":
            output.append("A")
        elif upper == "G":
            output.append("C")
        elif upper == "C":
            output.append("G")
        else:
            output.append("N")
    return "".join(output)


def cigar_contains_only_m(cigar: str) -> bool:
    """判断 CIGAR 字符串是否仅包含 M 操作（如 "100M"）"""
    return cigar.endswith("M") and cigar[:-1].isdigit()


def get_length_from_simple_cigar(cigar: str) -> int:
    """从简单 CIGAR（仅含 M）中提取长度"""
    return int(cigar[:-1])


def get_length_from_cigar(cigar: str) -> int:
    """
    从 CIGAR 字符串中计算重叠长度
    M/=/X/I/S/H/P 操作增加长度，D/N 操作减少长度
    """
    length = 0
    number: list[str] = []
    for char in cigar:
        if char.isdigit():
            number.append(char)
            continue
        if not number:
            continue
        count = int("".join(number))
        number.clear()
        if char in {"M", "=", "X", "I", "S", "H", "P"}:
            length += count
        elif char in {"D", "N"}:
            length -= count
    return length


def parse_int_tag(fields: list[str], prefix: str, fallback: int) -> int:
    """从 GFA 行的标签字段中解析整数值（如 "LN:i:1000" → 1000）"""
    for field in fields[3:]:
        if field.startswith(prefix):
            try:
                return int(field[len(prefix):])
            except ValueError:
                return fallback
    return fallback


def parse_float_tag(fields: list[str], prefix: str, fallback: float) -> float:
    """从 GFA 行的标签字段中解析浮点数值（如 "DP:f:30.5" → 30.5）"""
    for field in fields[3:]:
        if field.startswith(prefix):
            try:
                return float(field[len(prefix):])
            except ValueError:
                return fallback
    return fallback


def parse_string_tag(fields: list[str], *prefixes: str) -> str | None:
    """从 GFA 行的标签字段中解析字符串值（如 "ID:Z:edge1" → "edge1"）"""
    for prefix in prefixes:
        for field in fields:
            if field.startswith(prefix):
                return field[len(prefix):]
    return None


# ==================== 图视图层（用于渲染展示） ====================

@dataclass(eq=False)
class GraphNodeView:
    """
    图节点的视图对象
    从 DeBruijnNode 提取渲染所需的信息，用于画布显示
    """
    id: str                                              # 节点标识（不含方向标记）
    incident: list["GraphLinkView"] = field(default_factory=list)  # 与该节点关联的边
    line_points: list[Point2D] = field(default_factory=list)      # 节点路径的坐标点序列
    length: int = 1                                      # 序列长度
    depth: float = 0.0                                   # 测序深度
    width: float = 0.0                                   # 绘制宽度
    color: RgbaColor = field(default_factory=lambda: RgbaColor(128, 128, 128, 255))  # 节点颜色

    def first(self) -> Point2D:
        """路径起点"""
        return self.line_points[0]

    def second(self) -> Point2D:
        """路径第二个点（不足2个则返回起点）"""
        return self.line_points[1] if len(self.line_points) > 1 else self.line_points[0]

    def last(self) -> Point2D:
        """路径终点"""
        return self.line_points[-1]

    def second_last(self) -> Point2D:
        """路径倒数第二个点（不足2个则返回终点）"""
        return self.line_points[-2] if len(self.line_points) > 1 else self.line_points[-1]

    def centre_on_path(self) -> Point2D:
        """计算路径的中点（沿路径长度方向的中点，而非边界框中心）"""
        if not self.line_points:
            return Point2D(0.0, 0.0)
        if len(self.line_points) == 1:
            return self.line_points[0]
        total = 0.0
        for index in range(len(self.line_points) - 1):
            total += distance_points(self.line_points[index], self.line_points[index + 1])
        target = total / 2.0
        traversed = 0.0
        for index in range(len(self.line_points) - 1):
            a = self.line_points[index]
            b = self.line_points[index + 1]
            segment = distance_points(a, b)
            if traversed + segment >= target:
                fraction = 0.0 if segment < EPSILON else (target - traversed) / segment
                return Point2D(a.x + (b.x - a.x) * fraction, a.y + (b.y - a.y) * fraction)
            traversed += segment
        return self.line_points[len(self.line_points) // 2]


@dataclass(eq=False)
class GraphLinkView:
    """图边的视图对象，用于渲染展示"""
    from_node: GraphNodeView  # 起始节点视图
    to_node: GraphNodeView    # 终止节点视图
    from_orient: str          # 起始方向（+/-）
    to_orient: str            # 终止方向（+/-）
    read_count: int           # 重叠/read 数量
    name: str | None          # 边名称


# ==================== GFA 图视图（渲染管理层） ====================

class GfaGraph:
    """
    GFA 图视图类
    管理从 AssemblyGraph 到渲染层的转换，提供标签生成、视图同步等功能
    """
    EDGE_WIDTH = 1.5
    EDGE_LENGTH = 5.0

    def __init__(self, assembly_graph: AssemblyGraph) -> None:
        self.nodes: dict[str, GraphNodeView] = {}  # 渲染用节点视图（按 ID 索引）
        self.links: list[GraphLinkView] = []        # 渲染用边视图列表
        self.assembly_graph = assembly_graph         # 底层组装图

    @classmethod
    def parse(cls, input_path: Path, settings: BandageSettings, rng: random.Random | None = None) -> "GfaGraph":
        """从 GFA 文件解析并构建图视图"""
        return cls(AssemblyGraph.load_from_gfa(input_path, settings, rng))

    def prepare_styles(self) -> None:
        """准备渲染样式（计算深度、宽度、颜色）"""
        self.assembly_graph.prepare_for_rendering()

    def layout(self) -> None:
        """执行布局计算并同步视图数据"""
        self.assembly_graph.render_whole_graph()
        self.sync_from_assembly()

    def trim_depth(self, depth: float) -> str:
        """格式化深度值"""
        return format_trimmed_depth(depth)

    def node_label(self, node: GraphNodeView, options: BandageRenderOptions) -> str | None:
        """根据标签模式生成节点标签文本"""
        if options.node_label_mode == NodeLabelMode.NONE:
            return None
        if options.node_label_mode == NodeLabelMode.NAME:
            return node.id
        if options.node_label_mode == NodeLabelMode.DEPTH:
            return f"depth={self.trim_depth(node.depth)}"
        if options.node_label_mode == NodeLabelMode.LENGTH:
            return f"length={node.length}"
        if options.node_label_mode == NodeLabelMode.NAME_DEPTH:
            return f"{node.id}\ndepth={self.trim_depth(node.depth)}"
        if options.node_label_mode == NodeLabelMode.NAME_LENGTH:
            return f"{node.id}\nlength={node.length}"
        return f"{node.id}\ndepth={self.trim_depth(node.depth)}\nlength={node.length}"

    def edge_label(self, link: GraphLinkView, options: BandageRenderOptions) -> str | None:
        """根据标签模式生成边标签文本"""
        if options.edge_label_mode == EdgeLabelMode.NONE:
            return None
        if options.edge_label_mode == EdgeLabelMode.NAME:
            return link.name if link.name and link.name.strip() else None
        return f"{link.from_node.id}{link.from_orient} -> {link.to_node.id}{link.to_orient}"

    def sync_from_assembly(self) -> None:
        """
        从底层 AssemblyGraph 同步视图数据
        将 DeBruijnNode/Edge 转换为 GraphNodeView/GraphLinkView
        """
        self.nodes.clear()
        self.links.clear()
        view_nodes: dict[DeBruijnNode, GraphNodeView] = {}
        for node in self.assembly_graph.get_displayed_nodes():
            view_node = GraphNodeView(node.get_name_without_sign())
            view_node.length = node.length
            view_node.depth = node.depth
            view_node.width = node.width
            view_node.color = node.color
            view_node.line_points.extend(node.build_display_points())
            self.nodes[view_node.id] = view_node
            view_nodes[node] = view_node

        for edge in self.assembly_graph.get_displayed_edges():
            # 确定边的可见端节点（如果原始节点未绘制，使用其反向互补节点）
            starting_visible_node = edge.starting_node if edge.starting_node.drawn else edge.starting_node.reverse_complement
            ending_visible_node = edge.ending_node if edge.ending_node.drawn else edge.ending_node.reverse_complement
            if starting_visible_node is None or ending_visible_node is None:
                continue
            from_node = view_nodes.get(starting_visible_node)
            to_node = view_nodes.get(ending_visible_node)
            if from_node is None or to_node is None:
                continue
            link = GraphLinkView(
                from_node=from_node,
                to_node=to_node,
                from_orient="+" if edge.starting_node.drawn else "-",
                to_orient="+" if edge.ending_node.drawn else "-",
                read_count=edge.overlap,
                name=edge.label,
            )
            self.links.append(link)
            from_node.incident.append(link)
            to_node.incident.append(link)


# ==================== 矢量路径（SVG/PNG 渲染的几何路径） ====================

class VectorPath:
    """
    矢量路径，类似于 SVG 的 path 元素
    支持 MoveTo (M)、LineTo (L)、CurveTo (C)、Close (Z) 命令
    用于描述节点和边的绘制路径
    """
    def __init__(self) -> None:
        self.commands: list[tuple[str, tuple[float, ...]]] = []

    def move_to(self, x: float, y: float) -> None:
        """移动画笔到指定位置（M 命令）"""
        self.commands.append(("M", (x, y)))

    def line_to(self, x: float, y: float) -> None:
        """画直线到指定位置（L 命令）"""
        self.commands.append(("L", (x, y)))

    def curve_to(self, c1x: float, c1y: float, c2x: float, c2y: float, x: float, y: float) -> None:
        """画三次贝塞尔曲线到指定位置（C 命令），c1/c2 为控制点"""
        self.commands.append(("C", (c1x, c1y, c2x, c2y, x, y)))

    def close(self) -> None:
        """闭合路径（Z 命令）"""
        self.commands.append(("Z", ()))

    def sampled_points(self, flatness: float = 1.0) -> list[Point2D]:
        """
        将路径采样为离散点序列
        贝塞尔曲线根据其近似弧长自适应确定采样步数
        用于 PNG 渲染（Pillow 的 line 绘制需要点序列）和边界计算
        """
        points: list[Point2D] = []
        current = Point2D(0.0, 0.0)
        start = Point2D(0.0, 0.0)
        for command, values in self.commands:
            if command == "M":
                current = Point2D(values[0], values[1])
                start = current
                points.append(current)
            elif command == "L":
                current = Point2D(values[0], values[1])
                points.append(current)
            elif command == "C":
                c1 = Point2D(values[0], values[1])
                c2 = Point2D(values[2], values[3])
                end = Point2D(values[4], values[5])
                approx = distance_points(current, c1) + distance_points(c1, c2) + distance_points(c2, end)
                steps = max(12, min(96, int(max(approx / max(flatness, 0.25), 12))))
                for index in range(1, steps + 1):
                    points.append(cubic_point(current, c1, c2, end, index / steps))
                current = end
            elif command == "Z":
                points.append(start)
                current = start
        return points

    def svg_path(self, transform) -> str:
        """
        生成 SVG path 的 d 属性字符串
        transform 参数是一个坐标变换函数（世界坐标 → 屏幕坐标）
        """
        builder: list[str] = []
        for command, values in self.commands:
            if command == "M":
                point = transform(Point2D(values[0], values[1]))
                builder.append(f"M{format_number(point.x)} {format_number(point.y)}")
            elif command == "L":
                point = transform(Point2D(values[0], values[1]))
                builder.append(f"L{format_number(point.x)} {format_number(point.y)}")
            elif command == "C":
                c1 = transform(Point2D(values[0], values[1]))
                c2 = transform(Point2D(values[2], values[3]))
                end = transform(Point2D(values[4], values[5]))
                builder.append(
                    f"C{format_number(c1.x)} {format_number(c1.y)} "
                    f"{format_number(c2.x)} {format_number(c2.y)} "
                    f"{format_number(end.x)} {format_number(end.y)}"
                )
            elif command == "Z":
                builder.append("Z")
        return "".join(builder)


def cubic_point(p0: Point2D, p1: Point2D, p2: Point2D, p3: Point2D, t: float) -> Point2D:
    """计算三次贝塞尔曲线在参数 t 处的点"""
    one_minus_t = 1.0 - t
    x = one_minus_t ** 3 * p0.x + 3.0 * one_minus_t ** 2 * t * p1.x + 3.0 * one_minus_t * t ** 2 * p2.x + t ** 3 * p3.x
    y = one_minus_t ** 3 * p0.y + 3.0 * one_minus_t ** 2 * t * p1.y + 3.0 * one_minus_t * t ** 2 * p2.y + t ** 3 * p3.y
    return Point2D(x, y)


# ==================== 图形绘制项（节点和边的渲染路径） ====================

class BandageGraphicsItemNode:
    """节点的图形绘制项：将节点视图转换为矢量路径"""
    def __init__(self, node: GraphNodeView, settings: BandageSettings) -> None:
        self.node = node
        self.settings = settings
        self.path = VectorPath()
        self.remake_path()

    def get_first(self) -> Point2D:
        """路径起点"""
        return self.node.first()

    def get_second(self) -> Point2D:
        """路径第二个点"""
        return self.node.second()

    def get_last(self) -> Point2D:
        """路径终点"""
        return self.node.last()

    def get_second_last(self) -> Point2D:
        """路径倒数第二个点"""
        return self.node.second_last()

    def remake_path(self) -> None:
        """根据节点的坐标点序列构建矢量路径"""
        self.path = VectorPath()
        if not self.node.line_points:
            return
        first = self.node.line_points[0]
        self.path.move_to(first.x, first.y)
        for point in self.node.line_points[1:]:
            self.path.line_to(point.x, point.y)


class BandageGraphicsItemEdge:
    """
    边的图形绘制项：将边视图转换为带贝塞尔曲线的矢量路径
    边使用三次贝塞尔曲线连接节点端点，控制点沿节点方向延伸
    特殊情况：
    - 节点连接到自身反向互补：生成 S 形曲线路径
    - 节点连接到自身：生成环形曲线路径
    """
    def __init__(self, link: GraphLinkView, settings: BandageSettings, node_items: dict[GraphNodeView, BandageGraphicsItemNode]) -> None:
        self.link = link
        self.settings = settings
        self.node_items = node_items

    def path(self) -> VectorPath:
        """计算并返回边的矢量路径"""
        starting_location, before_starting_location, ending_location, after_ending_location = self.set_control_point_locations()
        edge_distance = distance_points(starting_location, ending_location)
        extension_length = min(self.settings.edge_length, edge_distance / 2.0)
        # 沿节点方向延伸控制点，使曲线与节点路径相切
        control_point1 = extend_line(before_starting_location, starting_location, extension_length)
        control_point2 = extend_line(after_ending_location, ending_location, extension_length)

        # 特殊情况：节点连接到其反向互补（如 1+ → 1-），在单链模式下
        if not self.settings.double_mode and self.link.from_node is self.link.to_node and self.link.from_orient != self.link.to_orient:
            reverse_extension = min(self.settings.edge_length / 2.0, edge_distance / 2.0)
            reverse_control1 = extend_line(before_starting_location, starting_location, reverse_extension)
            reverse_control2 = extend_line(after_ending_location, ending_location, reverse_extension)
            return make_special_path_connecting_node_to_reverse_complement(starting_location, ending_location, reverse_control1, reverse_control2)

        # 特殊情况：节点自环（如 1+ → 1+），生成环形路径
        if self.link.from_node is self.link.to_node:
            item = self.node_items.get(self.link.from_node)
            if item is not None and len(self.link.from_node.line_points) == 2:
                return make_special_path_connecting_node_to_self(starting_location, ending_location, control_point1, control_point2, self.settings.edge_length)

        # 一般情况：贝塞尔曲线路径
        path = VectorPath()
        path.move_to(starting_location.x, starting_location.y)
        path.curve_to(control_point1.x, control_point1.y, control_point2.x, control_point2.y, ending_location.x, ending_location.y)
        return path

    def set_control_point_locations(self) -> tuple[Point2D, Point2D, Point2D, Point2D]:
        """
        确定边的起止位置和控制点参考位置
        - 正链(+)方向：边从节点路径末端出发，到达目标节点路径首端
        - 负链(-)方向：边从节点路径首端出发，到达目标节点路径末端
        """
        starting_node = self.node_items.get(self.link.from_node)
        ending_node = self.node_items.get(self.link.to_node)
        if starting_node is None:
            raise RuntimeError(f"Missing node graphics item for {self.link.from_node.id}")
        if ending_node is None:
            raise RuntimeError(f"Missing node graphics item for {self.link.to_node.id}")
        starting_location = starting_node.get_last() if self.link.from_orient == "+" else starting_node.get_first()
        before_starting_location = starting_node.get_second_last() if self.link.from_orient == "+" else starting_node.get_second()
        ending_location = ending_node.get_first() if self.link.to_orient == "+" else ending_node.get_last()
        after_ending_location = ending_node.get_second() if self.link.to_orient == "+" else ending_node.get_second_last()
        return starting_location, before_starting_location, ending_location, after_ending_location


# ==================== 边路径辅助函数 ====================

def extend_line(start: Point2D, end: Point2D, extension_length: float) -> Point2D:
    """沿 start→end 方向延伸指定长度，返回延伸后的点"""
    dist = distance_points(start, end)
    if dist < 1e-6:
        return Point2D(end.x, end.y)
    ratio = extension_length / dist
    return Point2D(end.x + (end.x - start.x) * ratio, end.y + (end.y - start.y) * ratio)


def unit_normal(start: Point2D, end: Point2D) -> Point2D:
    """计算 start→end 方向的单位法向量（左手法则，逆时针旋转 90°）"""
    dx = end.x - start.x
    dy = end.y - start.y
    length = max(1e-6, math.hypot(dx, dy))
    return Point2D(-dy / length, dx / length)


def make_special_path_connecting_node_to_self(starting_location: Point2D, ending_location: Point2D, control_point1: Point2D, control_point2: Point2D, edge_length: float) -> VectorPath:
    """为节点自环边（如 1+ → 1+）创建环形曲线路径"""
    node_line_normal = unit_normal(starting_location, ending_location)
    perpendicular_shift = Point2D(node_line_normal.x * edge_length, node_line_normal.y * edge_length)
    node_mid_point = Point2D((starting_location.x + ending_location.x) / 2.0, (starting_location.y + ending_location.y) / 2.0)
    path = VectorPath()
    path.move_to(starting_location.x, starting_location.y)
    path.curve_to(
        control_point1.x,
        control_point1.y,
        control_point1.x + perpendicular_shift.x,
        control_point1.y + perpendicular_shift.y,
        node_mid_point.x + perpendicular_shift.x,
        node_mid_point.y + perpendicular_shift.y,
    )
    path.curve_to(
        control_point2.x + perpendicular_shift.x,
        control_point2.y + perpendicular_shift.y,
        control_point2.x,
        control_point2.y,
        ending_location.x,
        ending_location.y,
    )
    return path


def make_special_path_connecting_node_to_reverse_complement(starting_location: Point2D, ending_location: Point2D, control_point1: Point2D, control_point2: Point2D) -> VectorPath:
    """为节点连接到其反向互补（如 1+ → 1-）创建 S 形曲线路径"""
    start_to_control = Point2D(control_point1.x - starting_location.x, control_point1.y - starting_location.y)
    path_mid_point = Point2D(starting_location.x + start_to_control.x * 3.0, starting_location.y + start_to_control.y * 3.0)
    normal = unit_normal(control_point1, starting_location)
    normal_length = max(1e-6, distance_points(control_point1, starting_location))
    perpendicular_shift = Point2D(normal.x * normal_length * 1.5, normal.y * normal_length * 1.5)
    path = VectorPath()
    path.move_to(starting_location.x, starting_location.y)
    path.curve_to(
        control_point1.x,
        control_point1.y,
        path_mid_point.x + perpendicular_shift.x,
        path_mid_point.y + perpendicular_shift.y,
        path_mid_point.x,
        path_mid_point.y,
    )
    path.curve_to(
        path_mid_point.x - perpendicular_shift.x,
        path_mid_point.y - perpendicular_shift.y,
        control_point2.x,
        control_point2.y,
        ending_location.x,
        ending_location.y,
    )
    return path


# ==================== 渲染场景 ====================

@dataclass(frozen=True)
class NodeEntry:
    """渲染场景中的节点条目"""
    node: GraphNodeView         # 节点视图
    world_path: VectorPath      # 世界坐标系中的路径
    label_world_point: Point2D  # 标签位置（世界坐标）
    world_width: float          # 绘制宽度（世界坐标）


@dataclass(frozen=True)
class EdgeEntry:
    """渲染场景中的边条目"""
    link: GraphLinkView         # 边视图
    world_path: VectorPath      # 世界坐标系中的路径
    label_world_point: Point2D  # 标签位置（世界坐标）


@dataclass
class Bounds:
    """二维边界框"""
    min_x: float = math.inf
    min_y: float = math.inf
    max_x: float = -math.inf
    max_y: float = -math.inf

    def include(self, x: float, y: float) -> None:
        """扩展边界框以包含指定点"""
        self.min_x = min(self.min_x, x)
        self.min_y = min(self.min_y, y)
        self.max_x = max(self.max_x, x)
        self.max_y = max(self.max_y, y)

    def width(self) -> float:
        """边界框宽度"""
        return self.max_x - self.min_x

    def height(self) -> float:
        """边界框高度"""
        return self.max_y - self.min_y

    def pad(self, amount: float) -> None:
        """向四周扩展边界框指定距离"""
        if not math.isfinite(self.min_x):
            self.min_x = -amount
            self.min_y = -amount
            self.max_x = amount
            self.max_y = amount
        else:
            self.min_x -= amount
            self.min_y -= amount
            self.max_x += amount
            self.max_y += amount


class BandageRenderScene:
    """
    渲染场景：管理世界坐标到屏幕坐标的变换
    自动计算缩放比例和偏移量，使整图适配指定的画布尺寸
    """
    MARGIN = 96.0  # 画布边距（像素）

    def __init__(self, width: int, height: int, scale: float, offset_x: float, offset_y: float, node_items: dict[GraphNodeView, BandageGraphicsItemNode], nodes: list[NodeEntry], edges: list[EdgeEntry]) -> None:
        self.width = width       # 画布宽度
        self.height = height     # 画布高度
        self.scale = scale       # 缩放比例（世界坐标 → 屏幕坐标）
        self.offset_x = offset_x  # x 偏移量
        self.offset_y = offset_y  # y 偏移量
        self.node_items = node_items
        self.nodes = nodes
        self.edges = edges

    @classmethod
    def create(cls, graph: GfaGraph, width: int, height: int, settings: BandageSettings) -> "BandageRenderScene":
        """创建渲染场景，计算缩放和偏移以适配画布"""
        node_items: dict[GraphNodeView, BandageGraphicsItemNode] = {}
        for node in graph.nodes.values():
            node_items[node] = BandageGraphicsItemNode(node, settings)

        nodes: list[NodeEntry] = []
        for node_item in node_items.values():
            nodes.append(NodeEntry(node_item.node, node_item.path, node_item.node.centre_on_path(), node_item.node.width))

        edges: list[EdgeEntry] = []
        for link in graph.links:
            path = BandageGraphicsItemEdge(link, settings, node_items).path()
            edges.append(EdgeEntry(link, path, midpoint_of_path(path)))

        # 计算所有路径的边界框
        bounds = Bounds()
        for node in nodes:
            include_stroked_path(bounds, node.world_path, node.world_width)
        for edge in edges:
            include_stroked_path(bounds, edge.world_path, settings.edge_width)
        bounds.pad(24.0)

        # 计算缩放比例，使图适配画布（保持宽高比）
        world_width = max(bounds.width(), 1.0)
        world_height = max(bounds.height(), 1.0)
        available_width = max(1.0, width - cls.MARGIN * 2.0)
        available_height = max(1.0, height - cls.MARGIN * 2.0)
        scale = max(min(available_width / world_width, available_height / world_height), 0.1)
        # 计算偏移量，使图在画布中居中
        offset_x = (width - world_width * scale) / 2.0 - bounds.min_x * scale
        offset_y = (height - world_height * scale) / 2.0 - bounds.min_y * scale
        return cls(width, height, scale, offset_x, offset_y, node_items, nodes, edges)

    def to_screen_point(self, point: Point2D) -> Point2D:
        """世界坐标 → 屏幕坐标"""
        return Point2D(self.offset_x + point.x * self.scale, self.offset_y + point.y * self.scale)

    def scaled_width(self, world_width: float) -> float:
        """将世界坐标宽度转换为屏幕像素宽度（最小 0.75 像素）"""
        return max(world_width * self.scale, 0.75)


def include_stroked_path(bounds: Bounds, path: VectorPath, stroke_width: float) -> None:
    """将路径的描边区域纳入边界框计算"""
    radius = stroke_width / 2.0
    for point in path.sampled_points():
        bounds.include(point.x - radius, point.y - radius)
        bounds.include(point.x + radius, point.y + radius)


def midpoint_of_path(path: VectorPath) -> Point2D:
    """计算路径的中点（沿路径长度方向的中间位置）"""
    points = path.sampled_points()
    if not points:
        return Point2D(0.0, 0.0)
    if len(points) == 1:
        return points[0]
    total_length = 0.0
    for index in range(len(points) - 1):
        total_length += distance_points(points[index], points[index + 1])
    halfway = total_length / 2.0
    traversed = 0.0
    for index in range(len(points) - 1):
        a = points[index]
        b = points[index + 1]
        segment_length = distance_points(a, b)
        if traversed + segment_length >= halfway:
            fraction = 0.0 if segment_length < 1e-6 else (halfway - traversed) / segment_length
            return Point2D(a.x + (b.x - a.x) * fraction, a.y + (b.y - a.y) * fraction)
        traversed += segment_length
    return points[len(points) // 2]


# ==================== FM³ 力导向布局算法 ====================
# FM³ (Fast Multipole Multilevel Method) 是一种高效的多层级力导向图布局算法
# 核心思想：
#   1. 多层级粗化：将原图逐层压缩为更小的图（太阳系模型）
#   2. 自顶向下布局：从最粗的图开始布局，逐层细化
#   3. 力计算：使用弹簧吸引力 + 库仑斥力模型
#   4. 快速多极子方法（NMM）：用四叉树近似远距离斥力，降低 O(N²) 复杂度
#   5. 后处理：修复扭曲分裂、微调布局

@dataclass(eq=False)
class LevelNode:
    """
    多层级布局中的节点
    继承 LayoutNode 的所有力计算属性，并增加了多层级相关属性
    太阳系模型中的节点类型：
      - type 0: 未分类
      - type 1: 太阳节点（sun，被选为中心节点）
      - type 2: 行星节点（planet，与太阳直接相邻）
      - type 3: 行星-月球节点（planet-moon，同时有行星和月球属性）
      - type 4: 月球节点（moon，与行星相邻但不与太阳相邻）
    """
    id: int
    edges: list["LevelEdge"] = field(default_factory=list)
    original_node: LayoutNode | None = None
    lower_level_node: "LevelNode | None" = None
    higher_level_node: "LevelNode | None" = None
    x: float = 0.0
    y: float = 0.0
    width: float = 0.0
    height: float = 0.0
    force_x: float = 0.0
    force_y: float = 0.0
    last_move_x: float = 0.0
    last_move_y: float = 0.0
    attr_x: float = 0.0
    attr_y: float = 0.0
    rep_x: float = 0.0
    rep_y: float = 0.0
    mass: int = 0
    type: int = 0
    dedicated_sun_node: "LevelNode | None" = None
    dedicated_sun_distance: float = 0.0
    dedicated_pm_node: "LevelNode | None" = None
    lambda_list: list[float] = field(default_factory=list)
    neighbour_sun_node_list: list["LevelNode"] = field(default_factory=list)
    moon_nodes: list["LevelNode"] = field(default_factory=list)
    placed: bool = False
    angle1: float = 0.0
    angle2: float = TWO_PI

    def init_mult_values(self) -> None:
        """重置多层级相关的属性（在每一层粗化开始时调用）"""
        self.type = 0
        self.dedicated_sun_node = None
        self.dedicated_sun_distance = 0.0
        self.dedicated_pm_node = None
        self.lambda_list.clear()
        self.neighbour_sun_node_list.clear()
        self.moon_nodes.clear()
        self.placed = False
        self.angle1 = 0.0
        self.angle2 = TWO_PI


@dataclass(eq=False)
class LevelEdge:
    """多层级布局中的边"""
    source: LevelNode
    target: LevelNode
    length: float          # 期望长度
    moon_edge: bool = False  # 是否为月球边（连接月球节点和行星节点）
    extra_edge: bool = False # 是否为额外边

    def init_mult_values(self) -> None:
        """重置多层级相关的边属性"""
        self.moon_edge = False
        self.extra_edge = False


class LevelGraph:
    """多层级布局中的图结构"""
    def __init__(self) -> None:
        self.nodes: list[LevelNode] = []
        self.edges: list[LevelEdge] = []
        self.down_left_x = 0.0       # 四叉树包围盒左下角 x
        self.down_left_y = 0.0       # 四叉树包围盒左下角 y
        self.box_length = 0.0        # 四叉树包围盒边长
        self.average_ideal_edge_length = 0.0  # 平均理想边长

    def new_node(self) -> LevelNode:
        node = LevelNode(len(self.nodes))
        self.nodes.append(node)
        return node

    def new_edge(self, source: LevelNode, target: LevelNode, length: float) -> LevelEdge:
        edge = LevelEdge(source, target, length)
        self.edges.append(edge)
        source.edges.append(edge)
        target.edges.append(edge)
        return edge


@dataclass
class LayoutQuality:
    fixed_iterations: int
    fine_tuning_iterations: int


@dataclass
class FollowResult:
    finish: LevelNode
    path: list[LevelNode]
    steps: int


@dataclass
class PackingRectangle:
    width: float
    height: float
    old_down_left_x: float
    old_down_left_y: float
    component_index: int
    new_down_left_x: float = 0.0
    new_down_left_y: float = 0.0

    def area(self) -> float:
        return self.width * self.height


@dataclass
class PackingRow:
    max_height: float
    total_width: float
    row_index: int


@dataclass
class RectangleBounds:
    width: float
    height: float
    old_down_left_x: float
    old_down_left_y: float


@dataclass
class EdgeAccumulator:
    source: LevelNode
    target: LevelNode
    total_length: float = 0.0
    count: int = 0


@dataclass
class Component:
    nodes: list[LayoutNode]
    edges: list[LayoutEdge]
    min_x: float = 0.0
    min_y: float = 0.0
    max_x: float = 0.0
    max_y: float = 0.0
    width: float = 0.0
    height: float = 0.0

    def refresh_bounds(self, padding: float) -> None:
        self.min_x = math.inf
        self.min_y = math.inf
        self.max_x = -math.inf
        self.max_y = -math.inf
        for node in self.nodes:
            self.min_x = min(self.min_x, node.x)
            self.min_y = min(self.min_y, node.y)
            self.max_x = max(self.max_x, node.x)
            self.max_y = max(self.max_y, node.y)
        self.min_x -= padding / 2.0
        self.min_y -= padding / 2.0
        self.max_x += padding / 2.0
        self.max_y += padding / 2.0
        self.width = max(self.max_x - self.min_x, 0.0)
        self.height = max(self.max_y - self.min_y, 0.0)


def other(edge: LevelEdge, node: LevelNode) -> LevelNode:
    """获取边中除指定节点外的另一个端点"""
    return edge.target if edge.source is node else edge.source


class NodeSet:
    """
    节点集合，支持随机选择"星质量"最低的节点
    用于太阳系模型中选择太阳节点：优先选择邻居少的节点作为中心
    使用交换-删除策略实现 O(1) 的删除操作
    """
    def __init__(self, worker: "GraphLayoutWorker") -> None:
        self.worker = worker
        self.nodes: list[LevelNode] = []
        self.positions: dict[LevelNode, int] = {}
        self.star_mass: dict[LevelNode, int] = {}
        self.last_selectable_index = -1

    def init(self, graph: LevelGraph, weighted: bool) -> None:
        """初始化节点集合，weighted=True 时计算每个节点的星质量（自身+邻居质量之和）"""
        self.nodes.clear()
        self.positions.clear()
        self.star_mass.clear()
        for node in graph.nodes:
            self.positions[node] = len(self.nodes)
            self.nodes.append(node)
        self.last_selectable_index = len(self.nodes) - 1
        if weighted:
            for node in graph.nodes:
                mass = node.mass
                for edge in node.edges:
                    mass += other(edge, node).mass
                self.star_mass[node] = mass

    def is_empty(self) -> bool:
        """所有节点是否已被选完"""
        return self.last_selectable_index < 0

    def is_deleted(self, node: LevelNode) -> bool:
        """节点是否已被选走"""
        return self.positions[node] > self.last_selectable_index

    def delete_node(self, node: LevelNode) -> None:
        """删除节点（通过交换到末尾实现 O(1) 删除）"""
        delete_index = self.positions[node]
        last_selectable_node = self.nodes[self.last_selectable_index]
        self.nodes[self.last_selectable_index] = node
        self.nodes[delete_index] = last_selectable_node
        self.positions[node] = self.last_selectable_index
        self.positions[last_selectable_node] = delete_index
        self.last_selectable_index -= 1

    def get_random_node_with_lowest_star_mass(self, random_tries: int) -> LevelNode:
        """
        随机尝试 random_tries 次，返回星质量最低的节点
        这是太阳系模型中选择太阳节点的核心逻辑：
        随机采样若干候选，选最"孤立"的节点作为中心
        """
        selected_index = 0
        min_mass = 0
        random_node: LevelNode | None = None
        index = 1
        last_trie_index = self.last_selectable_index
        while index <= random_tries and last_trie_index >= 0:
            last_trie_node = self.nodes[last_trie_index]
            new_random_index = self.worker.random_int(0, last_trie_index)
            new_random_node = self.nodes[new_random_index]
            self.nodes[last_trie_index] = new_random_node
            self.nodes[new_random_index] = last_trie_node
            self.positions[new_random_node] = last_trie_index
            self.positions[last_trie_node] = new_random_index
            candidate = self.nodes[last_trie_index]
            if index == 1 or min_mass > self.star_mass[candidate]:
                selected_index = last_trie_index
                random_node = candidate
                min_mass = self.star_mass[candidate]
            index += 1
            last_trie_index -= 1
        assert random_node is not None
        last_selectable_node = self.nodes[self.last_selectable_index]
        self.nodes[self.last_selectable_index] = random_node
        self.nodes[selected_index] = last_selectable_node
        self.positions[random_node] = self.last_selectable_index
        self.positions[last_selectable_node] = selected_index
        self.last_selectable_index -= 1
        return random_node


# ==================== 四叉树（快速多极子方法的加速结构） ====================

class QuadCell:
    """
    四叉树单元格
    用于 NMM（快速多极子方法）中加速斥力计算
    核心思想：远距离节点群的斥力可以用质心近似，避免逐对计算
    - 叶节点存储实际的 LevelNode 列表
    - 非叶节点有 4 个子单元格
    - accumulate() 计算每个单元格的质心和总质量
    """
    def __init__(self, x: float, y: float, size: float, leaf_capacity: int) -> None:
        self.x = x
        self.y = y
        self.size = size
        self.leaf_capacity = leaf_capacity
        self.leaf_nodes: list[LevelNode] = []
        self.children: list["QuadCell | None"] | None = None
        self.mass = 0.0
        self.mass_x = 0.0
        self.mass_y = 0.0

    def insert(self, node: LevelNode) -> None:
        """插入节点到四叉树中；超出叶容量时自动细分"""
        if self.children is None:
            if len(self.leaf_nodes) < self.leaf_capacity or self.size < 1e-3:
                self.leaf_nodes.append(node)
                return
            nodes_to_fit = list(self.leaf_nodes)
            nodes_to_fit.append(node)
            self.shrink_to_small_cell(nodes_to_fit)
            self.subdivide()
            existing_nodes = list(self.leaf_nodes)
            self.leaf_nodes.clear()
            for existing in existing_nodes:
                self.insert_into_child(existing)
        self.insert_into_child(node)

    def accumulate(self) -> None:
        """递归计算每个单元格的质量和质心坐标"""
        if self.children is None:
            if not self.leaf_nodes:
                self.mass = 0.0
                self.mass_x = 0.0
                self.mass_y = 0.0
                return
            self.mass = float(len(self.leaf_nodes))
            self.mass_x = sum(node.x for node in self.leaf_nodes) / self.mass
            self.mass_y = sum(node.y for node in self.leaf_nodes) / self.mass
            return
        self.mass = 0.0
        self.mass_x = 0.0
        self.mass_y = 0.0
        for child in self.children:
            if child is None:
                continue
            child.accumulate()
            self.mass += child.mass
            self.mass_x += child.mass_x * child.mass
            self.mass_y += child.mass_y * child.mass
        if self.mass > 0.0:
            self.mass_x /= self.mass
            self.mass_y /= self.mass

    def collect_leaves(self, leaves: list["QuadCell"]) -> None:
        if self.children is None:
            if self.leaf_nodes:
                leaves.append(self)
            return
        for child in self.children:
            if child is not None:
                child.collect_leaves(leaves)

    def nodes(self) -> list[LevelNode]:
        return self.leaf_nodes

    def apply_repulsive_contribution(self, node: LevelNode, other_x: float, other_y: float, other_mass: float) -> None:
        """将一个质心对目标节点的斥力贡献施加到节点上（远距离近似计算）"""
        dx = node.x - other_x
        dy = node.y - other_y
        if abs(dx) < EPSILON and abs(dy) < EPSILON:
            return
        dist = math.hypot(dx, dy)
        if dist < EPSILON:
            return
        scalar = other_mass * (1.0 / dist) / dist
        node.rep_x += scalar * dx
        node.rep_y += scalar * dy

    def subdivide(self) -> None:
        """将单元格细分为 4 个子单元格"""
        self.children = [None, None, None, None]

    def insert_into_child(self, node: LevelNode) -> None:
        """将节点插入到对应的子单元格中（根据节点坐标判断象限）"""
        assert self.children is not None
        half = self.size / 2.0
        index = 0
        if node.x >= self.x + half:
            index += 1
        if node.y >= self.y + half:
            index += 2
        if self.children[index] is None:
            child_x = self.x if (index & 1) == 0 else self.x + half
            child_y = self.y if (index & 2) == 0 else self.y + half
            self.children[index] = QuadCell(child_x, child_y, half, self.leaf_capacity)
        self.children[index].insert(node)

    def shrink_to_small_cell(self, nodes_to_fit: list[LevelNode]) -> None:
        """将单元格缩小到刚好能容纳所有节点的最小象限（优化四叉树深度）"""
        while self.size > 1e-3:
            half = self.size / 2.0
            all_left = True
            all_right = True
            all_bottom = True
            all_top = True
            for node in nodes_to_fit:
                if node.x >= self.x + half:
                    all_left = False
                if node.x < self.x + half:
                    all_right = False
                if node.y >= self.y + half:
                    all_bottom = False
                if node.y < self.y + half:
                    all_top = False
            if all_left and all_bottom:
                self.size = half
            elif all_left and all_top:
                self.y += half
                self.size = half
            elif all_right and all_bottom:
                self.x += half
                self.size = half
            elif all_right and all_top:
                self.x += half
                self.y += half
                self.size = half
            else:
                break


class QuadTree:
    """四叉树，用于 NMM 斥力的快速近似计算"""
    def __init__(self, x: float, y: float, size: float, leaf_capacity: int) -> None:
        self.root = QuadCell(x, y, size, leaf_capacity)

    def insert(self, node: LevelNode) -> None:
        """插入节点到四叉树"""
        self.root.insert(node)

    def accumulate(self) -> None:
        """计算所有单元格的质心和质量"""
        self.root.accumulate()

    def collect_leaves(self) -> list[QuadCell]:
        """收集所有叶单元格"""
        leaves: list[QuadCell] = []
        self.root.collect_leaves(leaves)
        return leaves


# ==================== 快速多极子斥力计算 ====================

class NmmRepulsion:
    """
    快速多极子方法（NMM）斥力计算器
    节点数 >= NMM_MIN_NODE_NUMBER 时启用 NMM 近似，否则使用精确 O(N²) 计算
    NMM 斥力计算流程：
      1. 构建四叉树
      2. 计算每个叶单元格的质心
      3. 相邻叶单元格之间：逐对精确计算
      4. 远距离叶单元格之间：用质心近似计算
    """
    def __init__(self, worker: "GraphLayoutWorker") -> None:
        self.worker = worker
        self.using_nmm = False
        self.box_length = 0.0
        self.down_left_x = 0.0
        self.down_left_y = 0.0

    def make_initialisations(self, level_graph: LevelGraph) -> None:
        self.using_nmm = len(level_graph.nodes) >= GraphLayoutWorker.NMM_MIN_NODE_NUMBER
        self.box_length = level_graph.box_length
        self.down_left_x = level_graph.down_left_x
        self.down_left_y = level_graph.down_left_y

    def calculate_repulsive_forces(self, level_graph: LevelGraph) -> None:
        for node in level_graph.nodes:
            node.rep_x = 0.0
            node.rep_y = 0.0
        if self.using_nmm:
            self.calculate_repulsive_forces_by_nmm(level_graph)
        else:
            self.worker.calculate_exact_repulsive_forces(level_graph)

    def update_box_length_and_corner_coordinate(self, level_graph: LevelGraph) -> None:
        self.box_length = level_graph.box_length
        self.down_left_x = level_graph.down_left_x
        self.down_left_y = level_graph.down_left_y

    def deallocate_memory(self) -> None:
        self.using_nmm = False

    def calculate_repulsive_forces_by_nmm(self, level_graph: LevelGraph) -> None:
        tree = QuadTree(self.down_left_x, self.down_left_y, max(self.box_length, 1.0), GraphLayoutWorker.NMM_PARTICLES_IN_LEAVES)
        for node in level_graph.nodes:
            tree.insert(node)
        tree.accumulate()
        leaves = tree.collect_leaves()

        for leaf in leaves:
            self.add_direct_forces_inside_leaf(leaf)
        for index, act_leaf in enumerate(leaves):
            for other_leaf in leaves[index + 1:]:
                if self.bordering(act_leaf, other_leaf) or not self.well_separated(act_leaf, other_leaf):
                    self.add_direct_forces_between_leaves(act_leaf, other_leaf)
        for act_leaf in leaves:
            for other_leaf in leaves:
                if act_leaf is other_leaf or self.bordering(act_leaf, other_leaf) or not self.well_separated(act_leaf, other_leaf):
                    continue
                self.add_aggregated_leaf_force(other_leaf, act_leaf)

    def add_direct_forces_inside_leaf(self, leaf: QuadCell) -> None:
        nodes = leaf.nodes()
        for index, node in enumerate(nodes):
            for other_node in nodes[index + 1:]:
                self.worker.apply_repulsive_pair(node, other_node)

    def add_direct_forces_between_leaves(self, left: QuadCell, right: QuadCell) -> None:
        for left_node in left.nodes():
            for right_node in right.nodes():
                self.worker.apply_repulsive_pair(left_node, right_node)

    def add_aggregated_leaf_force(self, source_leaf: QuadCell, target_leaf: QuadCell) -> None:
        for node in target_leaf.nodes():
            source_leaf.apply_repulsive_contribution(node, source_leaf.mass_x, source_leaf.mass_y, source_leaf.mass)

    def bordering(self, node1: QuadCell, node2: QuadCell) -> bool:
        box_length1 = node1.size
        box_length2 = node2.size
        x1_min = node1.x
        x1_max = node1.x + box_length1
        y1_min = node1.y
        y1_max = node1.y + box_length1
        x2_min = node2.x
        x2_max = node2.x + box_length2
        y2_min = node2.y
        y2_max = node2.y + box_length2

        node2_contains1 = x2_min <= x1_min + EPSILON and x1_max <= x2_max + EPSILON and y2_min <= y1_min + EPSILON and y1_max <= y2_max + EPSILON
        node1_contains2 = x1_min <= x2_min + EPSILON and x2_max <= x1_max + EPSILON and y1_min <= y2_min + EPSILON and y2_max <= y1_max + EPSILON
        if node2_contains1 or node1_contains2:
            return False

        if box_length1 <= box_length2:
            if x1_min < x2_min:
                x1_min += box_length1
                x1_max += box_length1
            elif x1_max > x2_max:
                x1_min -= box_length1
                x1_max -= box_length1
            if y1_min < y2_min:
                y1_min += box_length1
                y1_max += box_length1
            elif y1_max > y2_max:
                y1_min -= box_length1
                y1_max -= box_length1
        else:
            if x2_min < x1_min:
                x2_min += box_length2
                x2_max += box_length2
            elif x2_max > x1_max:
                x2_min -= box_length2
                x2_max -= box_length2
            if y2_min < y1_min:
                y2_min += box_length2
                y2_max += box_length2
            elif y2_max > y1_max:
                y2_min -= box_length2
                y2_max -= box_length2

        shifted_node2_contains1 = x2_min <= x1_min + EPSILON and x1_max <= x2_max + EPSILON and y2_min <= y1_min + EPSILON and y1_max <= y2_max + EPSILON
        shifted_node1_contains2 = x1_min <= x2_min + EPSILON and x2_max <= x1_max + EPSILON and y1_min <= y2_min + EPSILON and y2_max <= y1_max + EPSILON
        return shifted_node2_contains1 or shifted_node1_contains2

    def well_separated(self, node1: QuadCell, node2: QuadCell) -> bool:
        box_length1 = node1.size
        box_length2 = node2.size
        if box_length1 <= box_length2:
            x1_min = node1.x
            x1_max = node1.x + box_length1
            y1_min = node1.y
            y1_max = node1.y + box_length1
            x2_min = node2.x - box_length2
            x2_max = node2.x + 2.0 * box_length2
            y2_min = node2.y - box_length2
            y2_max = node2.y + 2.0 * box_length2
        else:
            x1_min = node1.x - box_length1
            x1_max = node1.x + 2.0 * box_length1
            y1_min = node1.y - box_length1
            y1_max = node1.y + 2.0 * box_length1
            x2_min = node2.x
            x2_max = node2.x + box_length2
            y2_min = node2.y
            y2_max = node2.y + box_length2
        x_overlap = not (x1_max <= x2_min + EPSILON or x2_max <= x1_min + EPSILON)
        y_overlap = not (y1_max <= y2_min + EPSILON or y2_max <= y1_min + EPSILON)
        return not (x_overlap and y_overlap)


class Multilevel:
    """
    多层级粗化与展开模块
    实现 FM³ 算法的"太阳系"多层级图粗化策略：
    1. 将当前层的节点分组为"太阳系"（太阳+行星+月球）
    2. 将每个太阳系压缩为上一层的单个节点
    3. 在最粗层进行初始布局
    4. 逐层展开细化布局
    """
    def __init__(self, worker: "GraphLayoutWorker") -> None:
        self.worker = worker

    def create_multilevel_representations(self, level_zero: LevelGraph) -> list[LevelGraph]:
        """
        创建多层级图表示
        从第 0 层（原始图）开始，逐层粗化直到节点数 < MIN_GRAPH_SIZE 或边数不再显著减少
        """
        levels = [level_zero]
        bad_edge_counter = 0
        act_level = 0
        current = level_zero
        while len(current.nodes) > GraphLayoutWorker.MIN_GRAPH_SIZE and self.edge_numbers_of_all_levels_are_linear(levels, act_level, bad_edge_counter):
            next_level = LevelGraph()
            levels.append(next_level)
            self.init_multilevel_values(current)
            self.partition_galaxy_into_solar_systems(levels, act_level)  # 分组为太阳系
            self.collapse_solar_systems(levels, act_level)               # 压缩太阳系为上层节点
            act_level += 1
            current = levels[act_level]
            if len(current.edges) > 0.8 * len(levels[act_level - 1].edges):
                bad_edge_counter += 1
        return levels

    def find_initial_placement_for_level(self, level: int, levels: list[LevelGraph]) -> None:
        """为指定层级设置初始节点位置（从上层布局结果推算）"""
        pm_nodes: list[LevelNode] = []
        self.set_initial_positions_of_sun_nodes(level, levels)
        self.set_initial_positions_of_planet_and_moon_nodes(level, levels, pm_nodes)
        self.set_initial_positions_of_pm_nodes(level, levels, pm_nodes)

    def edge_numbers_of_all_levels_are_linear(self, levels: list[LevelGraph], act_level: int, bad_edge_counter: int) -> bool:
        if act_level == 0:
            return True
        current = levels[act_level]
        previous = levels[act_level - 1]
        if len(current.edges) <= 0.8 * len(previous.edges):
            return True
        return bad_edge_counter < 5

    def init_multilevel_values(self, level_graph: LevelGraph) -> None:
        for node in level_graph.nodes:
            node.init_mult_values()
        for edge in level_graph.edges:
            edge.init_mult_values()

    def partition_galaxy_into_solar_systems(self, levels: list[LevelGraph], act_level: int) -> None:
        self.create_suns_and_planets(levels, act_level)
        self.create_moon_nodes_and_pm_nodes(levels, act_level)

    def create_suns_and_planets(self, levels: list[LevelGraph], act_level: int) -> None:
        current = levels[act_level]
        next_level = levels[act_level + 1]
        node_set = NodeSet(self.worker)
        sun_nodes: list[LevelNode] = []

        for node in current.nodes:
            if act_level == 0:
                node.mass = 1
        node_set.init(current, True)

        while not node_set.is_empty():
            planet_nodes: list[LevelNode] = []
            sun_node = node_set.get_random_node_with_lowest_star_mass(GraphLayoutWorker.RANDOM_TRIES)
            sun_nodes.append(sun_node)

            new_node = next_level.new_node()
            sun_node.higher_level_node = new_node
            sun_node.type = 1
            sun_node.dedicated_sun_node = sun_node
            sun_node.dedicated_sun_distance = 0.0

            for edge in sun_node.edges:
                planet_node = other(edge, sun_node)
                planet_node.type = 2
                planet_node.dedicated_sun_node = sun_node
                planet_node.dedicated_sun_distance = edge.length
                planet_nodes.append(planet_node)

            for planet_node in planet_nodes:
                if not node_set.is_deleted(planet_node):
                    node_set.delete_node(planet_node)
            for planet_node in planet_nodes:
                for edge in planet_node.edges:
                    possible_moon_node = other(edge, planet_node)
                    if not node_set.is_deleted(possible_moon_node):
                        node_set.delete_node(possible_moon_node)

        for sun_node in sun_nodes:
            new_node = sun_node.higher_level_node
            if new_node is None:
                continue
            new_node.lower_level_node = sun_node
            new_node.x = sun_node.x
            new_node.y = sun_node.y
            new_node.width = sun_node.width
            new_node.height = sun_node.height
            new_node.mass = 0

    def create_moon_nodes_and_pm_nodes(self, levels: list[LevelGraph], act_level: int) -> None:
        current = levels[act_level]
        for node in current.nodes:
            if node.type != 0:
                continue
            nearest_neighbour_node: LevelNode | None = None
            dist_to_nearest_neighbour = 0.0
            moon_edge: LevelEdge | None = None
            first_adj_edge = True
            for edge in node.edges:
                neighbour_node = other(edge, node)
                if neighbour_node.type == 2 or neighbour_node.type == 3:
                    if first_adj_edge:
                        first_adj_edge = False
                        moon_edge = edge
                        dist_to_nearest_neighbour = edge.length
                        nearest_neighbour_node = neighbour_node
                    elif dist_to_nearest_neighbour > edge.length:
                        moon_edge = edge
                        dist_to_nearest_neighbour = edge.length
                        nearest_neighbour_node = neighbour_node
            if nearest_neighbour_node is None or moon_edge is None:
                continue
            moon_edge.moon_edge = True
            node.type = 4
            node.dedicated_sun_node = nearest_neighbour_node.dedicated_sun_node
            node.dedicated_sun_distance = dist_to_nearest_neighbour + nearest_neighbour_node.dedicated_sun_distance
            node.dedicated_pm_node = nearest_neighbour_node
            nearest_neighbour_node.type = 3
            nearest_neighbour_node.moon_nodes.append(node)

    def collapse_solar_systems(self, levels: list[LevelGraph], act_level: int) -> None:
        current = levels[act_level]
        next_level = levels[act_level + 1]
        for node in current.nodes:
            dedicated_sun = node.dedicated_sun_node
            if dedicated_sun is None or dedicated_sun.higher_level_node is None:
                continue
            dedicated_sun.higher_level_node.mass += 1

        merged: dict[tuple[int, int], EdgeAccumulator] = {}
        for edge in current.edges:
            s_node = edge.source
            t_node = edge.target
            s_sun_node = s_node.dedicated_sun_node
            t_sun_node = t_node.dedicated_sun_node
            if s_sun_node is None or t_sun_node is None or s_sun_node is t_sun_node:
                continue
            high_source = s_sun_node.higher_level_node
            high_target = t_sun_node.higher_level_node
            if high_source is None or high_target is None:
                continue
            new_length = s_node.dedicated_sun_distance + edge.length + t_node.dedicated_sun_distance
            key = edge_key(high_source.id, high_target.id)
            accumulator = merged.setdefault(key, EdgeAccumulator(high_source, high_target))
            accumulator.total_length += new_length
            accumulator.count += 1

            lambda_s = 0.0 if new_length == 0.0 else s_node.dedicated_sun_distance / new_length
            lambda_t = 0.0 if new_length == 0.0 else t_node.dedicated_sun_distance / new_length
            s_node.lambda_list.append(lambda_s)
            t_node.lambda_list.append(lambda_t)
            s_node.neighbour_sun_node_list.append(t_sun_node)
            t_node.neighbour_sun_node_list.append(s_sun_node)

        for accumulator in merged.values():
            next_level.new_edge(accumulator.source, accumulator.target, accumulator.total_length / accumulator.count)

    def set_initial_positions_of_sun_nodes(self, level: int, levels: list[LevelGraph]) -> None:
        higher = levels[level + 1]
        for higher_node in higher.nodes:
            lower_node = higher_node.lower_level_node
            if lower_node is None:
                continue
            lower_node.x = higher_node.x
            lower_node.y = higher_node.y
            lower_node.placed = True

    def set_initial_positions_of_planet_and_moon_nodes(self, level: int, levels: list[LevelGraph], pm_nodes: list[LevelNode]) -> None:
        current = levels[level]
        self.create_all_placement_sectors(level, levels)
        for node in current.nodes:
            if node.type == 3:
                pm_nodes.append(node)
            elif node.type == 2 or node.type == 4:
                positions: list[Point2D] = []
                dedicated_sun = node.dedicated_sun_node
                if dedicated_sun is None:
                    continue
                dedicated_sun_pos = Point2D(dedicated_sun.x, dedicated_sun.y)
                dedicated_sun_distance = node.dedicated_sun_distance

                for edge in node.edges:
                    adjacent = other(edge, node)
                    if node.dedicated_sun_node is adjacent.dedicated_sun_node and adjacent.type != 1 and adjacent.placed:
                        positions.append(self.calculate_position(dedicated_sun_pos, Point2D(adjacent.x, adjacent.y), dedicated_sun_distance, edge.length))

                if not node.lambda_list:
                    if not positions:
                        positions.append(self.create_random_pos(dedicated_sun_pos, node.dedicated_sun_distance, node.angle1, node.angle2))
                else:
                    for index, neighbour_sun_node in enumerate(node.neighbour_sun_node_list):
                        lambda_value = node.lambda_list[index % len(node.lambda_list)]
                        positions.append(self.get_waggled_inbetween_position(dedicated_sun_pos, Point2D(neighbour_sun_node.x, neighbour_sun_node.y), lambda_value))

                barycenter = self.get_barycenter_position(positions)
                node.x = barycenter.x
                node.y = barycenter.y
                node.placed = True

    def create_all_placement_sectors(self, level: int, levels: list[LevelGraph]) -> None:
        current = levels[level]
        higher = levels[level + 1]
        for higher_node in higher.nodes:
            adjacent_positions: list[Point2D] = []
            higher_pos = Point2D(higher_node.x, higher_node.y)
            for edge in higher_node.edges:
                if not edge.extra_edge:
                    adjacent = other(edge, higher_node)
                    adjacent_positions.append(Point2D(adjacent.x, adjacent.y))

            if not adjacent_positions:
                angle1 = 0.0
                angle2 = TWO_PI
            elif len(adjacent_positions) == 1:
                start_pos = adjacent_positions[0]
                x_parallel_pos = Point2D(higher_pos.x + 1.0, higher_pos.y)
                angle1 = angle_points(higher_pos, x_parallel_pos, start_pos)
                angle2 = angle1 + math.pi
            else:
                angle1 = 0.0
                angle2 = 0.0
                limit = min(len(adjacent_positions), 10)
                for index in range(limit):
                    start_pos = adjacent_positions[index]
                    x_parallel_pos = Point2D(higher_pos.x + 1.0, higher_pos.y)
                    act_angle1 = angle_points(higher_pos, x_parallel_pos, start_pos)
                    first_angle = True
                    min_next_angle = 0.0
                    for next_pos in adjacent_positions:
                        next_angle = angle_points(higher_pos, start_pos, next_pos)
                        if not same_point(start_pos, next_pos) and (first_angle or next_angle < min_next_angle):
                            min_next_angle = next_angle
                            first_angle = False
                    act_angle2 = act_angle1 + min_next_angle
                    if index == 0 or (act_angle2 - act_angle1) > (angle2 - angle1):
                        angle1 = act_angle1
                        angle2 = act_angle2
                if abs(angle1 - angle2) < EPSILON:
                    angle2 = angle1 + math.pi

            sun_node = higher_node.lower_level_node
            if sun_node is None:
                continue
            sun_node.angle1 = angle1
            sun_node.angle2 = angle2

        for node in current.nodes:
            dedicated_sun = node.dedicated_sun_node
            if dedicated_sun is not None:
                node.angle1 = dedicated_sun.angle1
                node.angle2 = dedicated_sun.angle2

    def set_initial_positions_of_pm_nodes(self, level: int, levels: list[LevelGraph], pm_nodes: list[LevelNode]) -> None:
        for node in pm_nodes:
            positions: list[Point2D] = []
            sun_node = node.dedicated_sun_node
            if sun_node is None:
                continue
            sun_pos = Point2D(sun_node.x, sun_node.y)
            sun_dist = node.dedicated_sun_distance

            for edge in node.edges:
                adjacent = other(edge, node)
                if not edge.moon_edge and node.dedicated_sun_node is adjacent.dedicated_sun_node and adjacent.type != 1 and adjacent.placed:
                    positions.append(self.calculate_position(sun_pos, Point2D(adjacent.x, adjacent.y), sun_dist, edge.length))

            for moon_node in node.moon_nodes:
                moon_pos = Point2D(moon_node.x, moon_node.y)
                moon_dist = moon_node.dedicated_sun_distance
                lambda_value = 0.0 if moon_dist == 0.0 else sun_dist / moon_dist
                positions.append(self.get_waggled_inbetween_position(sun_pos, moon_pos, lambda_value))

            if node.lambda_list:
                for index, neighbour_sun_node in enumerate(node.neighbour_sun_node_list):
                    lambda_value = node.lambda_list[index % len(node.lambda_list)]
                    positions.append(self.get_waggled_inbetween_position(sun_pos, Point2D(neighbour_sun_node.x, neighbour_sun_node.y), lambda_value))

            if not positions:
                positions.append(self.create_random_pos(sun_pos, sun_dist, node.angle1, node.angle2))

            barycenter = self.get_barycenter_position(positions)
            node.x = barycenter.x
            node.y = barycenter.y
            node.placed = True

    def create_random_pos(self, center: Point2D, radius: float, angle1: float, angle2: float) -> Point2D:
        rnd = (self.worker.random_int(1, 1_000_000_000) + 1.0) / 1_000_000_002.0
        rnd_angle = angle1 + (angle2 - angle1) * rnd
        return Point2D(center.x + math.cos(rnd_angle) * radius, center.y + math.sin(rnd_angle) * radius)

    def get_waggled_inbetween_position(self, source: Point2D, target: Point2D, lambda_value: float) -> Point2D:
        inbetween = Point2D(source.x + lambda_value * (target.x - source.x), source.y + lambda_value * (target.y - source.y))
        radius = GraphLayoutWorker.WAGGLE_FACTOR * distance_points(source, target)
        rnd = (self.worker.random_int(1, 1_000_000_000) + 1.0) / 1_000_000_002.0
        return self.create_random_pos(inbetween, radius * rnd, 0.0, TWO_PI)

    def get_barycenter_position(self, points: list[Point2D]) -> Point2D:
        if not points:
            return Point2D(0.0, 0.0)
        return Point2D(sum(point.x for point in points) / len(points), sum(point.y for point in points) / len(points))

    def calculate_position(self, p: Point2D, q: Point2D, dist_p: float, dist_q: float) -> Point2D:
        dist_pq = distance_points(p, q)
        if dist_pq < EPSILON:
            return self.create_random_pos(p, max(dist_p, dist_q), 0.0, TWO_PI)
        lambda_value = (dist_p + (dist_pq - dist_p - dist_q) / 2.0) / dist_pq
        return self.get_waggled_inbetween_position(p, q, lambda_value)


class GraphLayoutWorker:
    """
    图布局引擎：实现 FM³ 多层级力导向布局算法

    布局流程：
    1. 将图分解为连通分量
    2. 对每个分量独立布局：
       a. 构建多层级表示（粗化）
       b. 在最粗层进行随机初始布局
       c. 逐层展开并用力导向算法优化
    3. 旋转分量使包围盒面积最小
    4. 使用矩形装箱算法排列所有分量

    算法参数：
    - graph_layout_quality: 布局质量（0-4），值越高迭代次数越多
    - THRESHOLD: 力向量平均长度阈值，低于此值时停止迭代
    - FORCE_SCALING_FACTOR: 力的缩放因子
    - NMM_MIN_NODE_NUMBER: 启用快速多极子法的最小节点数
    """
    DEFAULT_ASPECT_RATIO = 1.333333
    THRESHOLD = 0.01
    FORCE_SCALING_FACTOR = 0.05
    FINE_TUNE_SCALAR = 0.2
    POST_SPRING_STRENGTH = 2.0
    STEPS_FOR_ROTATING_COMPONENTS = 50
    RESIZE_DRAWING = True
    RESIZING_SCALAR = 1.0
    MIN_NODE_SIZE = 10.0
    BOX_SCALING_FACTOR = 1.1
    MAX_ITER_FACTOR = 10
    MIN_GRAPH_SIZE = 50
    RANDOM_TRIES = 20
    WAGGLE_FACTOR = 0.05
    NMM_MIN_NODE_NUMBER = 175
    NMM_PARTICLES_IN_LEAVES = 25

    def __init__(self, graph: LayoutGraph, graph_layout_quality: int, linear_layout: bool, graph_layout_component_separation: float, aspect_ratio: float = DEFAULT_ASPECT_RATIO, rng: random.Random | None = None) -> None:
        self.graph = graph
        self.graph_layout_quality = graph_layout_quality
        self.linear_layout = linear_layout
        self.graph_layout_component_separation = graph_layout_component_separation
        self.aspect_ratio = aspect_ratio
        self.random = rng if rng is not None else random.Random(time.time_ns())
        self.nmm_repulsion = NmmRepulsion(self)

    def layout_graph(self) -> None:
        """布局整张图：分解连通分量 → 分别布局 → 旋转 → 装箱排列"""
        if not self.graph.nodes:
            return
        quality = self.quality(self.graph_layout_quality)
        components = self.connected_components()
        for component in components:
            self.layout_component(component, quality)
            self.rotate_component(component)
            component.refresh_bounds(self.graph_layout_component_separation)
        self.pack_components(components)

    def layout_component(self, component: Component, quality: LayoutQuality) -> None:
        """
        对单个连通分量执行 FM³ 布局
        1. 构建简化无环图
        2. 多层级粗化
        3. 自顶向下展开并力导向优化
        4. 导出最终坐标
        """
        if len(component.nodes) <= 1:
            component.refresh_bounds(self.graph_layout_component_separation)
            return
        level_zero = self.make_simple_loopfree(component)
        if len(level_zero.nodes) <= 1:
            self.export_positions(level_zero)
            component.refresh_bounds(self.graph_layout_component_separation)
            return
        multilevel = Multilevel(self)
        levels = multilevel.create_multilevel_representations(level_zero)
        max_level = len(levels) - 1
        for level in range(max_level, -1, -1):
            level_graph = levels[level]
            if level == max_level:
                self.create_initial_placement(level_graph, keep_positions=self.linear_layout)
            else:
                multilevel.find_initial_placement_for_level(level, levels)
                self.update_box_length_and_corner_coordinate(level_graph)
            self.call_force_calculation_step(level_graph, quality, level, max_level)
        self.export_positions(level_zero)
        component.refresh_bounds(0.0)

    def make_simple_loopfree(self, component: Component) -> LevelGraph:
        """构建简化无环图：移除自环边，合并重复边（取平均长度）"""
        reduced = LevelGraph()
        copies: dict[LayoutNode, LevelNode] = {}
        for original in component.nodes:
            copy = reduced.new_node()
            copy.original_node = original
            copy.x = original.x
            copy.y = original.y
            copies[original] = copy
        merged: dict[tuple[int, int], EdgeAccumulator] = {}
        for edge in component.edges:
            if edge.source is edge.target:
                continue
            source = copies[edge.source]
            target = copies[edge.target]
            key = edge_key(source.id, target.id)
            accumulator = merged.setdefault(key, EdgeAccumulator(source, target))
            accumulator.total_length += edge.desired_length
            accumulator.count += 1
        for accumulator in merged.values():
            reduced.new_edge(accumulator.source, accumulator.target, accumulator.total_length / accumulator.count)
        return reduced

    def export_positions(self, level_zero: LevelGraph) -> None:
        """将 LevelGraph 中的坐标导出到原始 LayoutNode"""
        for node in level_zero.nodes:
            if node.original_node is not None:
                node.original_node.x = node.x
                node.original_node.y = node.y

    def call_force_calculation_step(self, level_graph: LevelGraph, quality: LayoutQuality, act_level: int, max_level: int) -> None:
        """
        力计算的迭代主循环
        持续迭代直到：迭代次数达到上限 或 力向量平均长度低于阈值
        在第 0 层（最细节层）还会执行后处理步骤
        """
        if len(level_graph.nodes) <= 1:
            return
        iteration = 1
        max_mult_iter = self.get_max_mult_iter(act_level, max_level, quality.fixed_iterations, len(level_graph.nodes))
        average_force_length = self.THRESHOLD + 1.0
        self.set_average_ideal_edge_length(level_graph)
        self.nmm_repulsion.make_initialisations(level_graph)
        while iteration <= max_mult_iter and average_force_length >= self.THRESHOLD:
            self.calculate_forces(level_graph, iteration, 0, quality.fine_tuning_iterations)
            average_force_length = self.get_average_forcevector_length(level_graph)
            iteration += 1
        if act_level == 0:
            self.fix_twisted_splits(level_graph)
            self.call_post_processing_step(level_graph, quality.fine_tuning_iterations)
        self.nmm_repulsion.deallocate_memory()

    def call_post_processing_step(self, level_graph: LevelGraph, fine_tuning_iterations: int) -> None:
        """
        后处理步骤（仅在第 0 层执行）
        1. 10 次冷却迭代（cool_factor = 0.1）
        2. 缩放绘图使平均边长接近理想值
        3. fine_tuning_iterations 次微调迭代
        4. 再次缩放
        """
        for index in range(1, 11):
            self.calculate_forces(level_graph, index, 1, fine_tuning_iterations)
        if self.RESIZE_DRAWING:
            self.adapt_drawing_to_ideal_average_edge_length(level_graph)
            self.update_box_length_and_corner_coordinate(level_graph)
        for index in range(1, fine_tuning_iterations + 1):
            self.calculate_forces(level_graph, index, 2, fine_tuning_iterations)
        if self.RESIZE_DRAWING:
            self.adapt_drawing_to_ideal_average_edge_length(level_graph)
            self.update_box_length_and_corner_coordinate(level_graph)

    def calculate_forces(self, level_graph: LevelGraph, iteration: int, fine_tuning_step: int, fine_tuning_iterations: int) -> None:
        """单次力计算迭代：吸引力 → 斥力 → 合力 → 防振荡 → 移动节点"""
        self.calculate_attractive_forces(level_graph)
        self.nmm_repulsion.calculate_repulsive_forces(level_graph)
        self.add_attr_rep_forces(level_graph, iteration, fine_tuning_step, fine_tuning_iterations)
        self.prevent_oscillations(level_graph, iteration)
        self.move_nodes(level_graph)
        self.update_box_length_and_corner_coordinate(level_graph)

    def init_box_length_and_corner_coordinate(self, level_graph: LevelGraph) -> None:
        width_sum = 0.0
        height_sum = 0.0
        for node in level_graph.nodes:
            width_sum += max(node.width, self.MIN_NODE_SIZE)
            height_sum += max(node.height, self.MIN_NODE_SIZE)
        level_graph.box_length = math.ceil(max(width_sum, height_sum) * self.BOX_SCALING_FACTOR)
        level_graph.down_left_x = 0.0
        level_graph.down_left_y = 0.0

    def create_initial_placement(self, level_graph: LevelGraph, keep_positions: bool) -> None:
        self.init_box_length_and_corner_coordinate(level_graph)
        if keep_positions:
            self.update_box_length_and_corner_coordinate(level_graph)
            return
        for node in level_graph.nodes:
            node.x = 1.0 + self.random.random() * max(level_graph.box_length - 2.0, 1.0)
            node.y = 1.0 + self.random.random() * max(level_graph.box_length - 2.0, 1.0)
            node.last_move_x = 0.0
            node.last_move_y = 0.0
        self.update_box_length_and_corner_coordinate(level_graph)

    def calculate_attractive_forces(self, level_graph: LevelGraph) -> None:
        """计算所有边的弹簧吸引力（连接的节点相互吸引）"""
        for node in level_graph.nodes:
            node.attr_x = 0.0
            node.attr_y = 0.0
        for edge in level_graph.edges:
            source = edge.source
            target = edge.target
            dx = target.x - source.x
            dy = target.y - source.y
            dist = math.hypot(dx, dy)
            if dist < EPSILON:
                continue
            scalar = f_attr_scalar(dist, edge.length) / dist
            fx = scalar * dx
            fy = scalar * dy
            target.attr_x -= fx
            target.attr_y -= fy
            source.attr_x += fx
            source.attr_y += fy

    def calculate_exact_repulsive_forces(self, level_graph: LevelGraph) -> None:
        """精确计算所有节点对之间的斥力（O(N²)复杂度，仅用于小图）"""
        for index, source in enumerate(level_graph.nodes):
            for target in level_graph.nodes[index + 1:]:
                self.apply_repulsive_pair(source, target)

    def apply_repulsive_pair(self, source: LevelNode, target: LevelNode) -> None:
        """计算并施加一对节点之间的库仑斥力（反平方定律）"""
        dx = target.x - source.x
        dy = target.y - source.y
        if abs(dx) < EPSILON and abs(dy) < EPSILON:
            # 两个节点位置重叠时，用确定性伪随机方向推开
            angle = (((source.id * 1103515245 + target.id * 12345) & 1023) / 1024.0) * TWO_PI
            dx = math.cos(angle) * 1e-3
            dy = math.sin(angle) * 1e-3
        dist = math.hypot(dx, dy)
        if dist < EPSILON:
            return
        scalar = (1.0 / dist) / dist
        fx = scalar * dx
        fy = scalar * dy
        target.rep_x += fx
        target.rep_y += fy
        source.rep_x -= fx
        source.rep_y -= fy

    def add_attr_rep_forces(self, level_graph: LevelGraph, iteration: int, fine_tuning_step: int, fine_tuning_iterations: int) -> None:
        """
        合并吸引力和斥力，计算最终的节点合力
        考虑冷却因子（fine_tuning 阶段力会衰减）和最大移动半径限制
        """
        cool_factor = 1.0
        if fine_tuning_step == 1:
            cool_factor /= 10.0
        elif fine_tuning_step == 2:
            cool_factor = self.FINE_TUNE_SCALAR if iteration <= fine_tuning_iterations - 5 else self.FINE_TUNE_SCALAR / 10.0

        if fine_tuning_step <= 1:
            spring_strength = 1.0
            rep_strength = 1.0
        else:
            spring_strength = self.POST_SPRING_STRENGTH
            rep_strength = get_post_rep_force_strength(len(level_graph.nodes))

        ideal_square = level_graph.average_ideal_edge_length * level_graph.average_ideal_edge_length
        for node in level_graph.nodes:
            fx = (spring_strength * node.attr_x + rep_strength * node.rep_x) * ideal_square
            fy = (spring_strength * node.attr_y + rep_strength * node.rep_y) * ideal_square
            norm = math.hypot(fx, fy)
            if norm < EPSILON:
                node.force_x = 0.0
                node.force_y = 0.0
                continue
            scalar = min(norm * cool_factor * self.FORCE_SCALING_FACTOR, max_radius(level_graph.box_length, iteration)) / norm
            node.force_x = scalar * fx
            node.force_y = scalar * fy

    def prevent_oscillations(self, level_graph: LevelGraph, iteration: int) -> None:
        """
        防止节点振荡
        根据当前力方向与上次移动方向的夹角，限制节点移动距离
        夹角越小（方向一致），允许越大的移动
        夹角越大（方向反转），越强烈地限制移动
        """
        if iteration == 1:
            for node in level_graph.nodes:
                node.last_move_x = node.force_x
                node.last_move_y = node.force_y
            return

        pi1 = math.pi / 6.0
        pi2 = 2.0 * pi1
        pi3 = 3.0 * pi1
        pi4 = 4.0 * pi1
        pi5 = 5.0 * pi1
        pi7 = 7.0 * pi1
        pi8 = 8.0 * pi1
        pi9 = 9.0 * pi1
        pi10 = 10.0 * pi1
        pi11 = 11.0 * pi1

        for node in level_graph.nodes:
            norm_new = math.hypot(node.force_x, node.force_y)
            norm_old = math.hypot(node.last_move_x, node.last_move_y)
            if norm_new > EPSILON and norm_old > EPSILON:
                quotient = norm_old / norm_new
                angle_value = angle_vectors(node.last_move_x, node.last_move_y, node.force_x, node.force_y)
                factor = 1.0
                if (angle_value <= pi1 or angle_value >= pi11) and norm_new > norm_old * 2.0:
                    factor = quotient * 2.0
                elif pi1 <= angle_value <= pi2 and norm_new > norm_old * 1.5:
                    factor = quotient * 1.5
                elif pi2 <= angle_value <= pi3 and norm_new > norm_old:
                    factor = quotient
                elif pi3 <= angle_value <= pi4 and norm_new > norm_old * 0.66666666:
                    factor = quotient * 0.66666666
                elif pi4 <= angle_value <= pi5 and norm_new > norm_old * 0.5:
                    factor = quotient * 0.5
                elif pi5 <= angle_value <= pi7 and norm_new > norm_old * 0.33333333:
                    factor = quotient * 0.33333333
                elif pi7 <= angle_value <= pi8 and norm_new > norm_old * 0.5:
                    factor = quotient * 0.5
                elif pi8 <= angle_value <= pi9 and norm_new > norm_old * 0.66666666:
                    factor = quotient * 0.66666666
                elif pi9 <= angle_value <= pi10 and norm_new > norm_old:
                    factor = quotient
                elif pi10 <= angle_value <= pi11 and norm_new > norm_old * 1.5:
                    factor = quotient * 1.5
                node.force_x *= factor
                node.force_y *= factor
            node.last_move_x = node.force_x
            node.last_move_y = node.force_y

    def move_nodes(self, level_graph: LevelGraph) -> None:
        """根据计算好的力向量移动所有节点"""
        for node in level_graph.nodes:
            node.x += node.force_x
            node.y += node.force_y

    def update_box_length_and_corner_coordinate(self, level_graph: LevelGraph) -> None:
        """更新四叉树包围盒的范围（根据当前所有节点的坐标范围）"""
        first = level_graph.nodes[0]
        xmin = first.x
        xmax = first.x
        ymin = first.y
        ymax = first.y
        for node in level_graph.nodes:
            xmin = min(xmin, node.x)
            xmax = max(xmax, node.x)
            ymin = min(ymin, node.y)
            ymax = max(ymax, node.y)
        level_graph.down_left_x = math.floor(xmin - 1.0)
        level_graph.down_left_y = math.floor(ymin - 1.0)
        level_graph.box_length = math.ceil(max(ymax - ymin, xmax - xmin) * 1.01 + 2.0)
        if level_graph.box_length <= 2.0:
            level_graph.box_length = len(level_graph.nodes) * 20.0
            level_graph.down_left_x = math.floor(xmin) - level_graph.box_length / 2.0
            level_graph.down_left_y = math.floor(ymin) - level_graph.box_length / 2.0
        self.nmm_repulsion.update_box_length_and_corner_coordinate(level_graph)

    def set_average_ideal_edge_length(self, level_graph: LevelGraph) -> None:
        if not level_graph.edges:
            level_graph.average_ideal_edge_length = 50.0
            return
        level_graph.average_ideal_edge_length = sum(edge.length for edge in level_graph.edges) / len(level_graph.edges)

    def get_average_forcevector_length(self, level_graph: LevelGraph) -> float:
        return sum(math.hypot(node.force_x, node.force_y) for node in level_graph.nodes) / len(level_graph.nodes)

    def adapt_drawing_to_ideal_average_edge_length(self, level_graph: LevelGraph) -> None:
        sum_ideal = 0.0
        sum_real = 0.0
        for edge in level_graph.edges:
            sum_ideal += edge.length
            sum_real += math.hypot(edge.source.x - edge.target.x, edge.source.y - edge.target.y)
        area_scaling_factor = 1.0 if sum_real == 0.0 else sum_ideal / sum_real
        for node in level_graph.nodes:
            node.x = self.RESIZING_SCALAR * area_scaling_factor * node.x
            node.y = self.RESIZING_SCALAR * area_scaling_factor * node.y

    def fix_twisted_splits(self, level_graph: LevelGraph) -> None:
        """
        修复扭曲的分裂节点
        当一个节点有 3 个邻居，其中两条路径通向同一目标时，
        如果路径交叉则交换交叉点的坐标，使布局更清晰
        """
        for node in level_graph.nodes:
            adjacent_nodes = self.get_adjacent_nodes(node)
            if len(adjacent_nodes) != 3:
                continue

            direction1 = self.follow_nodes_until_branch(node, adjacent_nodes[0])
            direction2 = self.follow_nodes_until_branch(node, adjacent_nodes[1])
            direction3 = self.follow_nodes_until_branch(node, adjacent_nodes[2])
            path1: list[LevelNode] | None = None
            path2: list[LevelNode] | None = None
            if direction1.finish is direction2.finish and direction1.steps == direction2.steps and direction1.finish is not direction3.finish:
                path1 = direction1.path
                path2 = direction1.path
            elif direction1.finish is direction3.finish and direction1.steps == direction3.steps and direction1.finish is not direction2.finish:
                path1 = direction1.path
                path2 = direction3.path
            elif direction2.finish is direction3.finish and direction2.steps == direction3.steps and direction2.finish is not direction1.finish:
                path1 = direction2.path
                path2 = direction3.path
            if path1 is None or path2 is None or len(path1) <= 1 or len(path2) <= 1:
                continue
            for index in range(len(path1) - 1):
                path1_node1 = path1[index]
                path1_node2 = path1[index + 1]
                path2_node1 = path2[index]
                path2_node2 = path2[index + 1]
                if segment_intersects(
                    Point2D(path1_node1.x, path1_node1.y),
                    Point2D(path1_node2.x, path1_node2.y),
                    Point2D(path2_node1.x, path2_node1.y),
                    Point2D(path2_node2.x, path2_node2.y),
                ):
                    x = path1_node2.x
                    y = path1_node2.y
                    path1_node2.x = path2_node2.x
                    path1_node2.y = path2_node2.y
                    path2_node2.x = x
                    path2_node2.y = y

    def get_adjacent_nodes(self, node: LevelNode) -> list[LevelNode]:
        adjacent_nodes: list[LevelNode] = []
        for edge in node.edges:
            if edge.source is not node:
                adjacent_nodes.append(edge.source)
            if edge.target is not node:
                adjacent_nodes.append(edge.target)
        return adjacent_nodes

    def get_adjacent_nodes_excluding(self, node: LevelNode, excluded: LevelNode) -> list[LevelNode]:
        adjacent_nodes: list[LevelNode] = []
        for edge in node.edges:
            if edge.source is not node and edge.source is not excluded:
                adjacent_nodes.append(edge.source)
            if edge.target is not node and edge.source is not excluded:
                adjacent_nodes.append(edge.target)
        return adjacent_nodes

    def follow_nodes_until_branch(self, start: LevelNode, first: LevelNode) -> FollowResult:
        previous = start
        current = first
        steps = 0
        path: list[LevelNode] = []
        while True:
            adjacent_nodes = self.get_adjacent_nodes_excluding(current, previous)
            if len(adjacent_nodes) != 1:
                break
            previous = current
            current = adjacent_nodes[0]
            steps += 1
            path.append(previous)
        return FollowResult(current, path, steps)

    def rotate_component(self, component: Component) -> None:
        """
        旋转连通分量使其包围盒面积最小
        在 -45° ~ +45° 范围内尝试旋转，选择面积最小的角度
        如果最佳宽度 < 高度，额外旋转 90° 使宽 > 高
        """
        best_rectangle = self.calculate_bounding_rectangle(component)
        best_area = best_rectangle.width * best_rectangle.height
        best_coords = {node: Point2D(node.x, node.y) for node in component.nodes}
        old_coords = {node: Point2D(node.x, node.y) for node in component.nodes}
        pi4 = math.pi / 4.0

        for index in range(self.STEPS_FOR_ROTATING_COMPONENTS + 1):
            angle_value = (math.pi / 2.0) * (index / (self.STEPS_FOR_ROTATING_COMPONENTS + 1)) - pi4
            sin_value = math.sin(angle_value)
            cos_value = math.cos(angle_value)
            for node in component.nodes:
                point = old_coords[node]
                node.x = cos_value * point.x - sin_value * point.y
                node.y = sin_value * point.x + cos_value * point.y
            candidate = self.calculate_bounding_rectangle(component)
            area = candidate.width * candidate.height
            if area < best_area:
                best_area = area
                best_rectangle = candidate
                best_coords = {node: Point2D(node.x, node.y) for node in component.nodes}

        if best_rectangle.width / max(best_rectangle.height, EPSILON) < 1.0:
            for node, point in best_coords.items():
                node.x = -point.y
                node.y = point.x
        else:
            for node, point in best_coords.items():
                node.x = point.x
                node.y = point.y

    def calculate_bounding_rectangle(self, component: Component) -> RectangleBounds:
        xmin = xmax = ymin = ymax = 0.0
        first = True
        for node in component.nodes:
            act_x_min = act_x_max = node.x
            act_y_min = act_y_max = node.y
            if first:
                xmin = act_x_min
                xmax = act_x_max
                ymin = act_y_min
                ymax = act_y_max
                first = False
            else:
                xmin = min(xmin, act_x_min)
                xmax = max(xmax, act_x_max)
                ymin = min(ymin, act_y_min)
                ymax = max(ymax, act_y_max)
        xmin -= self.graph_layout_component_separation / 2.0
        xmax += self.graph_layout_component_separation / 2.0
        ymin -= self.graph_layout_component_separation / 2.0
        ymax += self.graph_layout_component_separation / 2.0
        return RectangleBounds(xmax - xmin, ymax - ymin, xmin, ymin)

    def pack_components(self, components: list[Component]) -> None:
        """使用矩形装箱算法将所有连通分量排列到画布上"""
        rectangles = [
            PackingRectangle(component.width, component.height, component.min_x, component.min_y, index)
            for index, component in enumerate(components)
        ]
        self.pack_rectangles_using_best_fit_strategy(rectangles)
        for rectangle in rectangles:
            component = components[rectangle.component_index]
            shift_x = rectangle.new_down_left_x - rectangle.old_down_left_x
            shift_y = rectangle.new_down_left_y - rectangle.old_down_left_y
            self.translate(component.nodes, shift_x, shift_y)
            component.refresh_bounds(0.0)

    def pack_rectangles_using_best_fit_strategy(self, rectangles: list[PackingRectangle]) -> None:
        """
        最佳适应矩形装箱算法
        1. 按面积从大到小排序
        2. 计算合适的换行宽度（使整体宽高比接近目标值）
        3. 逐行排列矩形
        """
        rectangles.sort(key=lambda rectangle: rectangle.area(), reverse=True)
        full_width = 0.0
        widest_rect = 0.0
        for rectangle in rectangles:
            full_width += rectangle.width
            widest_rect = max(widest_rect, rectangle.width)
        second_widest_rect = 0.0
        for rectangle in rectangles:
            if rectangle.width < widest_rect:
                second_widest_rect = max(second_widest_rect, rectangle.width)

        if second_widest_rect == 0.0:
            wrap_width = widest_rect
        elif widest_rect / second_widest_rect > 5.0:
            wrap_width = widest_rect
        else:
            wrap_width = full_width
            full_width_aspect_ratio = self.get_aspect_ratio(rectangles, full_width)
            best_agreement = self.get_aspect_ratio_agreement(self.aspect_ratio, full_width_aspect_ratio)
            if full_width_aspect_ratio > self.aspect_ratio:
                left = 0.0
                right = full_width
                while True:
                    mid = (left + right) / 2.0
                    mid_aspect_ratio = self.get_aspect_ratio(rectangles, mid)
                    if mid_aspect_ratio == self.aspect_ratio:
                        wrap_width = mid
                        break
                    elif mid_aspect_ratio > self.aspect_ratio:
                        right = mid
                    else:
                        left = mid
                    if wrap_width < widest_rect:
                        wrap_width = widest_rect
                        break
                    agreement = self.get_aspect_ratio_agreement(self.aspect_ratio, mid_aspect_ratio)
                    if agreement == best_agreement:
                        break
                    if agreement > best_agreement:
                        best_agreement = agreement
                        wrap_width = mid
                    if right - left < 1.0:
                        break

        rows: list[PackingRow] = []
        row_of_rectangle: list[int] = []
        width_of_current_row = 0.0
        for rectangle in rectangles:
            if not rows or width_of_current_row + rectangle.width > wrap_width or rectangle.width > wrap_width:
                rows.append(PackingRow(rectangle.height, rectangle.width, len(rows)))
                row_of_rectangle.append(len(rows) - 1)
                width_of_current_row = rectangle.width
            else:
                current_row = rows[-1]
                current_row.max_height = max(current_row.max_height, rectangle.height)
                current_row.total_width += rectangle.width
                row_of_rectangle.append(len(rows) - 1)
                width_of_current_row += rectangle.width
        self.export_new_rectangle_positions(rows, row_of_rectangle, rectangles)

    def export_new_rectangle_positions(self, rows: list[PackingRow], row_of_rectangle: list[int], rectangles: list[PackingRectangle]) -> None:
        row_y_min = [0.0] * len(rows)
        row_x_max = [0.0] * len(rows)
        for index in range(1, len(rows)):
            row_y_min[index] = row_y_min[index - 1] + rows[index - 1].max_height
        for index, rectangle in enumerate(rectangles):
            row = rows[row_of_rectangle[index]]
            new_x = row_x_max[row.row_index]
            row_x_max[row.row_index] += rectangle.width
            new_y = row_y_min[row.row_index] + (row.max_height - rectangle.height) / 2.0
            rectangle.new_down_left_x = new_x
            rectangle.new_down_left_y = new_y

    def get_aspect_ratio(self, rectangles: list[PackingRectangle], wrapping_width: float) -> float:
        width = height = row_width = row_height = 0.0
        for rectangle in rectangles:
            if row_width + rectangle.width <= wrapping_width:
                row_width += rectangle.width
                row_height = max(row_height, rectangle.height)
            else:
                width = max(width, row_width)
                height += row_height
                row_width = rectangle.width
                row_height = rectangle.height
        width = max(width, row_width)
        height += row_height
        return width / height if height > 0.0 else 1.0

    def get_aspect_ratio_agreement(self, ar1: float, ar2: float) -> float:
        if ar1 == 0.0 and ar2 == 0.0:
            return 1.0
        return min(ar1, ar2) / max(ar1, ar2)

    def connected_components(self) -> list[Component]:
        """使用 BFS 算法找出图中的所有连通分量"""
        components: list[Component] = []
        seen = [False] * len(self.graph.nodes)
        for node in self.graph.nodes:
            if seen[node.id]:
                continue
            component_nodes: list[LayoutNode] = []
            component_edges: list[LayoutEdge] = []
            queue: deque[LayoutNode] = deque([node])
            seen[node.id] = True
            while queue:
                current = queue.popleft()
                component_nodes.append(current)
                for edge in current.edges:
                    component_edges.append(edge)
                    other_node = edge.target if edge.source is current else edge.source
                    if not seen[other_node.id]:
                        seen[other_node.id] = True
                        queue.append(other_node)
            component = Component(component_nodes, self.deduplicate_edges(component_edges))
            component.refresh_bounds(0.0)
            components.append(component)
        return components

    def deduplicate_edges(self, edges: list[LayoutEdge]) -> list[LayoutEdge]:
        unique: dict[LayoutEdge, bool] = {}
        for edge in edges:
            unique[edge] = True
        return list(unique.keys())

    @staticmethod
    def translate(nodes: list[LayoutNode], shift_x: float, shift_y: float) -> None:
        for node in nodes:
            node.x += shift_x
            node.y += shift_y

    @staticmethod
    def quality(graph_layout_quality: int) -> LayoutQuality:
        """将布局质量等级（0-4）转换为具体的迭代次数配置"""
        if graph_layout_quality == 0:
            return LayoutQuality(3, 1)
        if graph_layout_quality == 1:
            return LayoutQuality(12, 8)
        if graph_layout_quality == 3:
            return LayoutQuality(60, 20)
        if graph_layout_quality == 4:
            return LayoutQuality(120, 20)
        return LayoutQuality(30, 20)  # 默认等级 2

    @staticmethod
    def get_max_mult_iter(act_level: int, max_level: int, fixed_iterations: int, node_count: int) -> int:
        """
        计算当前层级的最大迭代次数
        层级越高（图越粗），迭代次数越少
        小图（<=500 节点）至少 100 次迭代
        """
        if max_level == 0:
            iterations = fixed_iterations + (GraphLayoutWorker.MAX_ITER_FACTOR - 1) * fixed_iterations
        else:
            iterations = fixed_iterations + int((act_level / max(max_level, 1)) * (GraphLayoutWorker.MAX_ITER_FACTOR - 1) * fixed_iterations)
        return 100 if node_count <= 500 and iterations < 100 else iterations

    def random_int(self, min_inclusive: int, max_inclusive: int) -> int:
        if max_inclusive <= min_inclusive:
            return min_inclusive
        return self.random.randint(min_inclusive, max_inclusive)


# ==================== 渲染输出函数 ====================

def load_font(size: int, bold: bool) -> Any:
    """
    加载字体文件用于 PNG 渲染
    优先使用 Windows 系统字体（Arial），回退到默认字体
    """
    from PIL import ImageFont

    font_candidates: list[Path] = []
    windows_fonts = Path("C:/Windows/Fonts")
    if bold:
        font_candidates.extend([windows_fonts / "arialbd.ttf", windows_fonts / "Arial Bold.ttf", windows_fonts / "ARIALBD.TTF"])
    else:
        font_candidates.extend([windows_fonts / "arial.ttf", windows_fonts / "Arial.ttf", windows_fonts / "ARIAL.TTF"])
    font_candidates.extend([windows_fonts / "calibri.ttf", windows_fonts / "DejaVuSans.ttf"])
    for path in font_candidates:
        if path.exists():
            try:
                return ImageFont.truetype(str(path), size=size)
            except OSError:
                continue
    try:
        return ImageFont.truetype("arialbd.ttf" if bold else "arial.ttf", size=size)
    except OSError:
        return ImageFont.load_default()


def screen_points(scene: BandageRenderScene, path: VectorPath) -> list[tuple[float, float]]:
    """将路径的世界坐标点转换为屏幕坐标点（用于 Pillow 绘制）"""
    return [(scene.to_screen_point(point).x, scene.to_screen_point(point).y) for point in path.sampled_points()]


def draw_centered_label(draw: Any, text: str | None, point: Point2D, font: Any, fill: RgbaColor) -> None:
    """
    在 PNG 图像上绘制居中标签文字
    使用白色光晕（描边）提高文字在彩色背景上的可读性
    """
    if text is None or not text.strip():
        return
    lines = text.splitlines()
    if not lines:
        return
    bboxes = [draw.textbbox((0, 0), line, font=font) for line in lines]
    line_heights = [bottom - top for left, top, right, bottom in bboxes]
    line_spacing = max(2, int(round(font.size * 0.2))) if hasattr(font, "size") else 2
    total_height = sum(line_heights) + line_spacing * (len(lines) - 1)
    current_y = point.y - total_height / 2.0
    for line, bbox, line_height in zip(lines, bboxes, line_heights):
        left, top, right, bottom = bbox
        x = point.x - (right - left) / 2.0
        y = current_y
        for dx in range(-1, 2):
            for dy in range(-1, 2):
                if dx == 0 and dy == 0:
                    continue
                draw.text((x + dx, y + dy), line, font=font, fill=LABEL_HALO_COLOUR.to_pillow())
        draw.text((x, y), line, font=font, fill=fill.to_pillow())
        current_y += line_height + line_spacing


def render_png(graph: GfaGraph, output: Path, width: int, height: int, settings: BandageSettings, options: BandageRenderOptions) -> None:
    """
    将图渲染为 PNG 图片
    绘制顺序：边 → 节点 → 边标签 → 节点标签（标签在最上层）
    """
    try:
        from PIL import Image, ImageDraw
    except ImportError as exc:
        raise SystemExit("PNG rendering requires Pillow. Install it with 'pip install pillow'.") from exc

    if output.parent:
        output.parent.mkdir(parents=True, exist_ok=True)
    scene = BandageRenderScene.create(graph, width, height, settings)
    image = Image.new("RGBA", (width, height), settings.background_colour.to_pillow())
    draw = ImageDraw.Draw(image, "RGBA")

    # 绘制边
    edge_width = max(1, int(round(scene.scaled_width(settings.edge_width))))
    for edge in scene.edges:
        points = screen_points(scene, edge.world_path)
        if len(points) >= 2:
            draw.line(points, fill=settings.edge_colour.to_pillow(), width=edge_width, joint="curve")

    # 绘制节点
    for node in scene.nodes:
        points = screen_points(scene, node.world_path)
        node_width = max(1, int(round(scene.scaled_width(node.world_width))))
        if len(points) >= 2:
            draw.line(points, fill=node.node.color.to_pillow(), width=node_width, joint="curve")

    # 绘制标签
    edge_font = load_font(options.edge_font_size, bold=False)
    node_font = load_font(options.node_font_size, bold=True)
    for edge in scene.edges:
        draw_centered_label(draw, graph.edge_label(edge.link, options), scene.to_screen_point(edge.label_world_point), edge_font, options.font_color)
    for node in scene.nodes:
        draw_centered_label(draw, graph.node_label(node.node, options), scene.to_screen_point(node.label_world_point), node_font, options.font_color)

    image.save(output, format="PNG")


def write_svg_text(handle, text: str | None, point: Point2D, font_size: int, weight: str, fill: RgbaColor) -> None:
    """向 SVG 文件写入一个居中文本元素，带白色描边光晕"""
    if text is None or not text.strip():
        return
    lines = text.splitlines()
    line_spacing = font_size * 1.2
    start_y = point.y - (line_spacing * (len(lines) - 1)) / 2.0
    handle.write(
        f'    <text x="{format_number(point.x)}" y="{format_number(start_y)}" text-anchor="middle" dominant-baseline="middle" '
        f'font-family="Arial, sans-serif" font-size="{font_size}" font-weight="{weight}" fill="{fill.hex()}"'
    )
    if fill.a < 255:
        handle.write(f' fill-opacity="{format_number(fill.opacity())}"')
    handle.write(f' stroke="{LABEL_HALO_COLOUR.hex()}"')
    if LABEL_HALO_COLOUR.a < 255:
        handle.write(f' stroke-opacity="{format_number(LABEL_HALO_COLOUR.opacity())}"')
    handle.write(' stroke-width="3" paint-order="stroke">')
    for index, line in enumerate(lines):
        dy = "0" if index == 0 else format_number(line_spacing)
        handle.write(f'<tspan x="{format_number(point.x)}" dy="{dy}">{escape(line)}</tspan>')
    handle.write("</text>\n")


def render_svg(graph: GfaGraph, output: Path, width: int, height: int, settings: BandageSettings, options: BandageRenderOptions) -> None:
    """
    将图渲染为 SVG 文件
    使用 SVG 的 path 和 text 元素，支持无损缩放
    """
    if output.parent:
        output.parent.mkdir(parents=True, exist_ok=True)
    scene = BandageRenderScene.create(graph, width, height, settings)
    with output.open("w", encoding="utf-8") as handle:
        handle.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        handle.write(f'<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="{width}" height="{height}" viewBox="0 0 {width} {height}">\n')
        # 背景
        handle.write(f'  <rect x="0" y="0" width="100%" height="100%" fill="{settings.background_colour.hex()}"')
        if settings.background_colour.a < 255:
            handle.write(f' fill-opacity="{format_number(settings.background_colour.opacity())}"')
        handle.write('/>\n')
        # 边路径
        handle.write("  <g>\n")
        for edge in scene.edges:
            handle.write(f'    <path d="{edge.world_path.svg_path(scene.to_screen_point)}" fill="none" stroke="{settings.edge_colour.hex()}"')
            if settings.edge_colour.a < 255:
                handle.write(f' stroke-opacity="{format_number(settings.edge_colour.opacity())}"')
            handle.write(f' stroke-width="{format_number(scene.scaled_width(settings.edge_width))}" stroke-linecap="round" stroke-linejoin="round"/>\n')
        handle.write("  </g>\n")
        # 节点路径
        handle.write("  <g>\n")
        for node in scene.nodes:
            handle.write(f'    <path d="{node.world_path.svg_path(scene.to_screen_point)}" fill="none" stroke="{node.node.color.hex()}"')
            if node.node.color.a < 255:
                handle.write(f' stroke-opacity="{format_number(node.node.color.opacity())}"')
            handle.write(f' stroke-width="{format_number(scene.scaled_width(node.world_width))}" stroke-linecap="butt" stroke-linejoin="round"/>\n')
        handle.write("  </g>\n")
        # 标签
        handle.write("  <g>\n")
        for edge in scene.edges:
            write_svg_text(handle, graph.edge_label(edge.link, options), scene.to_screen_point(edge.label_world_point), options.edge_font_size, "400", options.font_color)
        for node in scene.nodes:
            write_svg_text(handle, graph.node_label(node.node, options), scene.to_screen_point(node.label_world_point), options.node_font_size, "700", options.font_color)
        handle.write("  </g>\n")
        handle.write("</svg>\n")


# ==================== 命令行入口 ====================

def build_parser() -> argparse.ArgumentParser:
    """构建命令行参数解析器"""
    parser = argparse.ArgumentParser(description="Pure Python Bandage-style GFA renderer. Parses GFA, computes layout, and renders PNG/SVG without Java.")
    parser.add_argument("input_gfa", help="Path to the input GFA file")
    parser.add_argument("output_image", help="Output image path (.png or .svg)")
    parser.add_argument("width", nargs="?", type=int, help="Optional output width")
    parser.add_argument("height", nargs="?", type=int, help="Optional output height")
    parser.add_argument("--node-label", choices=[mode.value for mode in NodeLabelMode], help="Node label mode")
    parser.add_argument("--edge-label", choices=[mode.value for mode in EdgeLabelMode], help="Edge label mode")
    parser.add_argument("--font-size", type=int, help="Set node and edge label font size")
    parser.add_argument("--node-font-size", type=int, help="Set node label font size")
    parser.add_argument("--edge-font-size", type=int, help="Set edge label font size")
    parser.add_argument("--font-color", help="Set label color, e.g. black or #000000")
    parser.add_argument("--seed", type=int, help="Optional random seed for stable layout/colors")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    """验证命令行参数的合法性"""
    if (args.width is None) != (args.height is None):
        raise SystemExit("Width and height must be provided together, or both omitted.")
    for value, name in ((args.font_size, "font-size"), (args.node_font_size, "node-font-size"), (args.edge_font_size, "edge-font-size")):
        if value is not None and value <= 0:
            raise SystemExit(f"{name} must be greater than zero.")


def build_render_options(args: argparse.Namespace) -> BandageRenderOptions:
    """从命令行参数构建渲染选项"""
    options = BandageRenderOptions()
    node_label_mode = options.node_label_mode if args.node_label is None else BandageRenderOptions.parse_node_label_mode(args.node_label)
    edge_label_mode = options.edge_label_mode if args.edge_label is None else BandageRenderOptions.parse_edge_label_mode(args.edge_label)
    node_font_size = options.node_font_size
    edge_font_size = options.edge_font_size
    if args.font_size is not None:
        node_font_size = args.font_size
        edge_font_size = args.font_size
    if args.node_font_size is not None:
        node_font_size = args.node_font_size
    if args.edge_font_size is not None:
        edge_font_size = args.edge_font_size
    font_color = options.font_color if args.font_color is None else parse_color(args.font_color)
    return BandageRenderOptions(node_label_mode, edge_label_mode, node_font_size, edge_font_size, font_color)


def run(input_gfa: Path, output_image: Path, width: int, height: int, options: BandageRenderOptions, seed: int | None = None) -> None:
    """
    主渲染函数：解析 GFA → 布局计算 → 输出图片
    根据输出文件后缀选择 PNG 或 SVG 格式
    """
    settings = BandageSettings()
    rng = random.Random(seed if seed is not None else time.time_ns())
    graph = GfaGraph.parse(input_gfa, settings, rng)
    graph.prepare_styles()
    graph.layout()
    suffix = output_image.suffix.lower()
    if suffix == ".svg":
        render_svg(graph, output_image, width, height, settings, options)
    elif suffix == ".png":
        render_png(graph, output_image, width, height, settings, options)
    else:
        raise SystemExit(f"Unsupported output format: {output_image}")


# ==================== HiMT 集成接口 ====================

def _resolve_himt_output_path(args: argparse.Namespace) -> Path:
    """根据 HiMT 参数解析输出文件路径（输出目录 + 文件名 + 格式后缀）"""
    image_format = args.image_format.lower()
    if args.output_name:
        output_name = Path(args.output_name)
        if output_name.is_absolute():
            raise ValueError("--output_name must be a file name or a relative path inside the output directory.")
        if output_name.suffix and output_name.suffix.lower() != f".{image_format}":
            raise ValueError("The suffix of --output_name must match --image_format.")
        if not output_name.suffix:
            output_name = output_name.with_suffix(f".{image_format}")
    else:
        output_name = Path(args.input_file).stem
        output_name = Path(f"{output_name}.{image_format}")
    return Path(args.output_dir) / output_name


def _build_render_args(args: argparse.Namespace) -> argparse.Namespace:
    """从 HiMT 参数构建标准渲染参数"""
    return argparse.Namespace(
        width=args.width,
        height=args.height,
        node_label=args.node_label,
        edge_label=args.edge_label,
        font_size=args.font_size,
        node_font_size=args.node_font_size,
        edge_font_size=args.edge_font_size,
        font_color=args.font_color,
        seed=args.seed,
    )


def view(args: argparse.Namespace) -> None:
    """
    HiMT 集成的可视化入口函数
    接收 HiMT 的参数命名空间，执行 GFA 可视化渲染
    """
    input_path = Path(args.input_file).expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"Input GFA file not found: {args.input_file}")
    if args.width <= 0 or args.height <= 0:
        raise ValueError("Image width and height must be greater than zero.")

    render_args = _build_render_args(args)
    validate_args(render_args)
    output_path = _resolve_himt_output_path(args).expanduser().resolve()
    options = build_render_options(render_args)

    print(f"({datetime.datetime.now()}) rendering GFA view")
    print(f"({datetime.datetime.now()}) input gfa: {input_path}")
    print(f"({datetime.datetime.now()}) output image: {output_path}")

    run(input_path, output_path, args.width, args.height, options, args.seed)

    print(f"({datetime.datetime.now()}) view finished")


def main(argv: list[str] | None = None) -> int:
    """独立运行时的命令行入口"""
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(args)
    width = 2200 if args.width is None else args.width
    height = 1600 if args.height is None else args.height
    options = build_render_options(args)
    run(Path(args.input_gfa).expanduser().resolve(), Path(args.output_image).expanduser().resolve(), width, height, options, args.seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
