#!/usr/bin/env python3
"""
Geometry and export helpers for kirigami laser-pattern generators.

This module provides the shared Point/Polyline geometry types, SVG/DXF/PNG
writers, and the verified ULB-style tilted square/triangular sheet functions
used by the active generators in this folder.

The geometry follows the triangular formulas in the supplementary material of:
Rafsanjani and Pasini, "Bistable auxetic mechanical metamaterials inspired by
ancient geometric motifs", Extreme Mechanics Letters 9 (2016) 291-296.

For triangular designs:
    tilted:   a/l = cos(theta) - sqrt(3) sin(theta)
    circular: a/l = (-1 + sqrt(12 (R/l)^2 - 3)) / 2
    parallel: a/l = 1 - 2 sqrt(3) (w/l)

The script generates cut centerlines for triangular building blocks in a
rhombic triangular grid. For the tilted motif, each block is built from the
three mutually intersecting tilted slit lines shown in the paper figure; the
central rotating triangle is formed by their intersections. By default, the
output uses uninterrupted solid cut centerlines for laser software. If you want
to leave physical ligaments like the paper's finite hinge thickness t, set
`--hinge` to a positive value in millimeters.

Example:
    python triangular_bistable_auxetic.py --motif tilted --a-over-l 0.5 \
        --l 20 --hinge 0 --rows 2 --cols 2 --svg tri_bam.svg --dxf tri_bam.dxf
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple


SQRT3 = math.sqrt(3.0)
EPS = 1.0e-9


@dataclass(frozen=True)
class Point:
    x: float
    y: float

    def __add__(self, other: "Point") -> "Point":
        return Point(self.x + other.x, self.y + other.y)

    def __sub__(self, other: "Point") -> "Point":
        return Point(self.x - other.x, self.y - other.y)

    def __mul__(self, scalar: float) -> "Point":
        return Point(self.x * scalar, self.y * scalar)

    __rmul__ = __mul__

    def dot(self, other: "Point") -> float:
        return self.x * other.x + self.y * other.y

    def cross(self, other: "Point") -> float:
        return self.x * other.y - self.y * other.x

    def norm(self) -> float:
        return math.hypot(self.x, self.y)

    def unit(self) -> "Point":
        n = self.norm()
        if n < EPS:
            raise ValueError("Cannot normalize a zero-length vector")
        return Point(self.x / n, self.y / n)

    def left_normal(self) -> "Point":
        return Point(-self.y, self.x)

    def rotate(self, angle: float) -> "Point":
        c = math.cos(angle)
        s = math.sin(angle)
        return Point(c * self.x - s * self.y, s * self.x + c * self.y)


@dataclass(frozen=True)
class Polyline:
    points: Tuple[Point, ...]
    layer: str = "CUT"

    def length(self) -> float:
        return sum(distance(a, b) for a, b in zip(self.points[:-1], self.points[1:]))


@dataclass(frozen=True)
class PatternParameters:
    motif: str
    l: float
    a_over_l: float
    theta: float | None = None
    radius_over_l: float | None = None
    width_over_l: float | None = None

    @property
    def a(self) -> float:
        return self.a_over_l * self.l

    @property
    def full_expansion_strain(self) -> float:
        if self.motif == "tilted":
            if self.theta is None:
                raise ValueError("tilted motif requires theta")
            return 2.0 * self.a_over_l * math.sin(self.theta + math.pi / 6.0)
        return self.a_over_l


def distance(a: Point, b: Point) -> float:
    return (b - a).norm()


def line_intersection(p1: Point, d1: Point, p2: Point, d2: Point) -> Point:
    den = d1.cross(d2)
    if abs(den) < EPS:
        raise ValueError("Parallel lines do not intersect")
    t = (p2 - p1).cross(d2) / den
    return p1 + d1 * t


def point_line_parameter(origin: Point, direction: Point, point: Point) -> float:
    return (point - origin).dot(direction) / direction.dot(direction)


def point_almost_equal(a: Point, b: Point, tol: float = 1.0e-7) -> bool:
    return distance(a, b) <= tol


def segment_line_intersection(a: Point, b: Point, origin: Point, direction: Point) -> Point | None:
    edge = b - a
    den = direction.cross(edge)
    if abs(den) < EPS:
        return None
    t = (a - origin).cross(edge) / den
    u = (a - origin).cross(direction) / den
    if -1.0e-8 <= u <= 1.0 + 1.0e-8:
        return origin + direction * t
    return None


def clipped_line_in_triangle(
    origin: Point,
    direction: Point,
    triangle: Sequence[Point],
) -> Tuple[Point, Point]:
    return clipped_line_in_polygon(origin, direction, ensure_ccw(triangle))


def clipped_line_in_polygon(
    origin: Point,
    direction: Point,
    polygon: Sequence[Point],
) -> Tuple[Point, Point]:
    verts = tuple(polygon)
    hits: List[Point] = []
    for i in range(len(verts)):
        hit = segment_line_intersection(verts[i], verts[(i + 1) % len(verts)], origin, direction)
        if hit is not None and not any(point_almost_equal(hit, p) for p in hits):
            hits.append(hit)

    if len(hits) < 2:
        raise ValueError("Cut line does not cross the polygon")

    hits.sort(key=lambda p: point_line_parameter(origin, direction, p))
    return hits[0], hits[-1]


def split_line_with_gaps(
    p0: Point,
    p1: Point,
    gap_centers: Sequence[Point],
    endpoint_gap: float,
    intersection_gap: float,
) -> List[Polyline]:
    d = (p1 - p0).unit()
    total = distance(p0, p1)
    intervals = [(0.0, endpoint_gap), (total - endpoint_gap, total)]

    for center in gap_centers:
        s = (center - p0).dot(d)
        intervals.append((s - 0.5 * intersection_gap, s + 0.5 * intersection_gap))

    clipped = []
    for start, end in intervals:
        start = max(0.0, min(total, start))
        end = max(0.0, min(total, end))
        if end > start + EPS:
            clipped.append((start, end))
    clipped.sort()

    merged: List[Tuple[float, float]] = []
    for start, end in clipped:
        if not merged or start > merged[-1][1] + EPS:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))

    cut_intervals = []
    cursor = 0.0
    for start, end in merged:
        if start > cursor + EPS:
            cut_intervals.append((cursor, start))
        cursor = max(cursor, end)
    if total > cursor + EPS:
        cut_intervals.append((cursor, total))

    paths = []
    for start, end in cut_intervals:
        if end - start > EPS:
            paths.append(Polyline((p0 + d * start, p0 + d * end)))
    return paths


def angle_of(v: Point) -> float:
    return math.atan2(v.y, v.x)


def positive_angle_delta(a0: float, a1: float) -> float:
    delta = a1 - a0
    while delta < 0.0:
        delta += 2.0 * math.pi
    while delta >= 2.0 * math.pi:
        delta -= 2.0 * math.pi
    return delta


def shortest_angle_delta(a0: float, a1: float) -> float:
    delta = (a1 - a0 + math.pi) % (2.0 * math.pi) - math.pi
    return delta


def point_on_circle(center: Point, radius: float, angle: float) -> Point:
    return Point(center.x + radius * math.cos(angle), center.y + radius * math.sin(angle))


def arc_polyline(
    center: Point,
    radius: float,
    start: Point,
    end: Point,
    ccw: bool | None = True,
    segments: int = 24,
) -> Polyline:
    a0 = angle_of(start - center)
    a1 = angle_of(end - center)
    if ccw is None:
        delta = shortest_angle_delta(a0, a1)
    elif ccw:
        delta = positive_angle_delta(a0, a1)
    else:
        delta = -positive_angle_delta(a1, a0)

    n = max(2, int(math.ceil(abs(delta) / (math.pi / max(3, segments)))))
    pts = []
    for i in range(n + 1):
        t = i / float(n)
        pts.append(point_on_circle(center, radius, a0 + delta * t))
    return Polyline(tuple(pts))


def trim_polyline(polyline: Polyline, trim_start: float, trim_end: float) -> Polyline:
    points = list(polyline.points)
    total = polyline.length()
    if total <= trim_start + trim_end + EPS:
        raise ValueError(
            "Cut path is shorter than requested hinge trimming. "
            "Reduce hinge or increase a/l."
        )

    def point_at_distance(points_: Sequence[Point], s: float) -> Point:
        if s <= 0.0:
            return points_[0]
        remaining = s
        for p0, p1 in zip(points_[:-1], points_[1:]):
            seg = distance(p0, p1)
            if remaining <= seg + EPS:
                t = max(0.0, min(1.0, remaining / seg))
                return p0 + (p1 - p0) * t
            remaining -= seg
        return points_[-1]

    start = point_at_distance(points, trim_start)
    end = point_at_distance(points, total - trim_end)

    kept = [start]
    walked = 0.0
    for p0, p1 in zip(points[:-1], points[1:]):
        seg = distance(p0, p1)
        next_walked = walked + seg
        if next_walked > trim_start + EPS and walked < total - trim_end - EPS:
            if walked >= trim_start - EPS and next_walked <= total - trim_end + EPS:
                kept.append(p1)
        walked = next_walked
    if distance(kept[-1], end) > EPS:
        kept.append(end)
    if len(kept) < 2:
        kept = [start, end]
    return Polyline(tuple(kept), layer=polyline.layer)


def theta_from_a_over_l(a_over_l: float) -> float:
    if not 0.0 <= a_over_l <= 1.0:
        raise ValueError("For triangular tilted motifs, a/l must be in [0, 1]")
    return math.acos(a_over_l / 2.0) - math.pi / 3.0


def radius_over_l_from_a_over_l(a_over_l: float) -> float:
    if not 0.0 <= a_over_l <= 1.0:
        raise ValueError("For triangular circular motifs, a/l must be in [0, 1]")
    return math.sqrt(((2.0 * a_over_l + 1.0) ** 2 + 3.0) / 12.0)


def width_over_l_from_a_over_l(a_over_l: float) -> float:
    if not 0.0 <= a_over_l <= 1.0:
        raise ValueError("For triangular parallel motifs, a/l must be in [0, 1]")
    return (1.0 - a_over_l) / (2.0 * SQRT3)


def square_tilted_a_over_l(theta: float) -> float:
    if not 0.0 <= theta <= math.pi / 4.0 + EPS:
        raise ValueError("Square tilted theta must be in [0, 45] degrees")
    return math.cos(theta) - math.sin(theta)


def square_tilted_theta_from_a_over_l(a_over_l: float) -> float:
    if not 0.0 <= a_over_l <= 1.0:
        raise ValueError("For square tilted motifs, a/l must be in [0, 1]")
    return math.acos(a_over_l / math.sqrt(2.0)) - math.pi / 4.0


def derive_parameters(
    motif: str,
    l: float,
    a_over_l: float | None,
    theta_deg: float | None,
    radius_over_l: float | None,
    width_over_l: float | None,
) -> PatternParameters:
    motif = motif.lower()
    if motif not in {"tilted", "circular", "parallel"}:
        raise ValueError("motif must be one of: tilted, circular, parallel")
    if l <= 0.0:
        raise ValueError("l must be positive")

    if motif == "tilted":
        if theta_deg is not None:
            theta = math.radians(theta_deg)
            if not 0.0 <= theta <= math.pi / 6.0 + EPS:
                raise ValueError("Triangular tilted theta must be in [0, 30] degrees")
            actual_a = math.cos(theta) - SQRT3 * math.sin(theta)
        else:
            actual_a = 0.5 if a_over_l is None else a_over_l
            theta = theta_from_a_over_l(actual_a)
        return PatternParameters(motif, l, actual_a, theta=theta)

    if motif == "circular":
        if radius_over_l is not None:
            lo = SQRT3 / 3.0
            hi = SQRT3 / 2.0
            if not lo - EPS <= radius_over_l <= hi + EPS:
                raise ValueError("Triangular circular R/l must be in [sqrt(3)/3, sqrt(3)/2]")
            actual_a = (-1.0 + math.sqrt(12.0 * radius_over_l * radius_over_l - 3.0)) / 2.0
        else:
            actual_a = 0.5 if a_over_l is None else a_over_l
            radius_over_l = radius_over_l_from_a_over_l(actual_a)
        return PatternParameters(motif, l, actual_a, radius_over_l=radius_over_l)

    if width_over_l is not None:
        lo = 0.0
        hi = SQRT3 / 6.0
        if not lo - EPS <= width_over_l <= hi + EPS:
            raise ValueError("Triangular parallel w/l must be in [0, sqrt(3)/6]")
        actual_a = 1.0 - 2.0 * SQRT3 * width_over_l
    else:
        actual_a = 0.5 if a_over_l is None else a_over_l
        width_over_l = width_over_l_from_a_over_l(actual_a)
    return PatternParameters(motif, l, actual_a, width_over_l=width_over_l)


def ensure_ccw(triangle: Sequence[Point]) -> Tuple[Point, Point, Point]:
    if len(triangle) != 3:
        raise ValueError("A triangular building block requires exactly three vertices")
    a, b, c = triangle
    if (b - a).cross(c - a) < 0.0:
        return a, c, b
    return a, b, c


def tilted_paths(
    triangle: Sequence[Point],
    params: PatternParameters,
    hinge: float,
) -> List[Polyline]:
    verts = ensure_ccw(triangle)
    theta = params.theta
    if theta is None:
        raise ValueError("tilted motif requires theta")

    lines = []
    for i in range(3):
        p = verts[i]
        q = verts[(i + 1) % 3]
        d = (q - p).unit().rotate(theta)
        lines.append((p, d))

    # Vertex i of the rotating triangle is the intersection of line i-1 and i.
    inner = []
    for i in range(3):
        p0, d0 = lines[(i - 1) % 3]
        p1, d1 = lines[i]
        inner.append(line_intersection(p0, d0, p1, d1))

    paths: List[Polyline] = []
    for i in range(3):
        p, d = lines[i]
        start, end = clipped_line_in_triangle(p, d, verts)
        points_on_line = []
        for vertex in inner:
            s = point_line_parameter(start, end - start, vertex)
            if -1.0e-7 <= s <= 1.0 + 1.0e-7:
                points_on_line.append(vertex)
        paths.extend(
            split_line_with_gaps(
                start,
                end,
                points_on_line,
                endpoint_gap=0.5 * hinge,
                intersection_gap=hinge,
            )
        )
    return paths


def tilted_selected_paths(
    triangle: Sequence[Point],
    params: PatternParameters,
    ligament: float,
) -> List[Polyline]:
    verts = ensure_ccw(triangle)
    theta = params.theta
    if theta is None:
        raise ValueError("tilted motif requires theta")

    lines = []
    for i in range(3):
        p = verts[i]
        q = verts[(i + 1) % 3]
        d = (q - p).unit().rotate(theta)
        lines.append((p, d))

    inner = []
    for i in range(3):
        p0, d0 = lines[(i - 1) % 3]
        p1, d1 = lines[i]
        inner.append(line_intersection(p0, d0, p1, d1))

    paths: List[Polyline] = []
    for p, d in lines:
        start, end = clipped_line_in_polygon(p, d, verts)
        line_vec = end - start
        denom = line_vec.dot(line_vec)
        points_on_line = []
        for vertex in inner:
            s = (vertex - start).dot(line_vec) / denom
            if -1.0e-7 <= s <= 1.0 + 1.0e-7 and abs((vertex - start).cross(line_vec)) <= 1.0e-6:
                points_on_line.append(vertex)
        if len(points_on_line) != 2:
            continue

        # The laser motif shown in the paper keeps the long slit that starts
        # from the construction-line vertex and stops near the farther
        # rotating-unit corner. It is trimmed at both ends by t/2.
        endpoint = start if distance(start, p) <= distance(end, p) else end
        far_corner = max(points_on_line, key=lambda q: distance(q, endpoint))
        raw = Polyline((endpoint, far_corner))
        paths.append(trim_polyline(raw, 0.5 * ligament, 0.5 * ligament))
    return paths


def parallel_paths(
    triangle: Sequence[Point],
    params: PatternParameters,
    hinge: float,
) -> List[Polyline]:
    verts = ensure_ccw(triangle)
    if params.width_over_l is None:
        raise ValueError("parallel motif requires width_over_l")
    offset = params.width_over_l * params.l

    lines = []
    for i in range(3):
        p = verts[i]
        q = verts[(i + 1) % 3]
        d = (q - p).unit()
        inward = d.left_normal()
        lines.append((p + inward * offset, d))

    inner = []
    for i in range(3):
        p0, d0 = lines[(i - 1) % 3]
        p1, d1 = lines[i]
        inner.append(line_intersection(p0, d0, p1, d1))

    paths = []
    for i in range(3):
        raw = Polyline((inner[i], inner[(i + 1) % 3]))
        paths.append(trim_polyline(raw, hinge, hinge))
    return paths


def circular_paths(
    triangle: Sequence[Point],
    params: PatternParameters,
    hinge: float,
    arc_segments: int,
) -> List[Polyline]:
    verts = ensure_ccw(triangle)
    if params.radius_over_l is None:
        raise ValueError("circular motif requires radius_over_l")
    radius = params.radius_over_l * params.l

    side_points = []
    for i in range(3):
        p = verts[i]
        q = verts[(i + 1) % 3]
        side = q - p
        side_len = side.norm()
        if radius < side_len / 2.0:
            raise ValueError("Circle radius is too small for adjacent side intersections")
        midpoint = p + side * 0.5
        inward = side.unit().left_normal()
        height = math.sqrt(max(0.0, radius * radius - (side_len / 2.0) ** 2))
        side_points.append(midpoint + inward * height)

    paths = []
    for i in range(3):
        center = verts[i]
        start = side_points[i]
        end = side_points[(i - 1) % 3]
        raw = arc_polyline(center, radius, start, end, ccw=None, segments=arc_segments)
        paths.append(trim_polyline(raw, hinge, hinge))
    return paths


def square_tilted_cell_paths(
    origin: Point,
    cell_size: float,
    theta: float,
    ligament: float,
    mirror: bool = False,
) -> List[Polyline]:
    if cell_size <= 0.0:
        raise ValueError("cell_size must be positive")
    if ligament < 0.0:
        raise ValueError("t/ligament must be non-negative")

    verts = (
        origin,
        origin + Point(cell_size, 0.0),
        origin + Point(cell_size, cell_size),
        origin + Point(0.0, cell_size),
    )
    signed_theta = -theta if mirror else theta

    lines = []
    for i in range(4):
        p = verts[i]
        q = verts[(i + 1) % 4]
        d = (q - p).unit().rotate(signed_theta)
        lines.append((p, d))

    inner = []
    for i in range(4):
        p0, d0 = lines[(i - 1) % 4]
        p1, d1 = lines[i]
        inner.append(line_intersection(p0, d0, p1, d1))

    paths: List[Polyline] = []
    for p, d in lines:
        try:
            start, end = clipped_line_in_polygon(p, d, verts)
        except ValueError:
            continue
        denom = (end - start).dot(end - start)
        points_on_line = []
        for vertex in inner:
            s = (vertex - start).dot(end - start) / denom
            if -1.0e-7 <= s <= 1.0 + 1.0e-7:
                points_on_line.append(vertex)
        paths.extend(
            split_line_with_gaps(
                start,
                end,
                points_on_line,
                endpoint_gap=0.5 * ligament,
                intersection_gap=ligament,
            )
        )
    return paths


def clip_segment_to_rectangle(p0: Point, p1: Point, width: float, height: float) -> Tuple[Point, Point] | None:
    dx = p1.x - p0.x
    dy = p1.y - p0.y
    u0 = 0.0
    u1 = 1.0

    for p, q in ((-dx, p0.x), (dx, width - p0.x), (-dy, p0.y), (dy, height - p0.y)):
        if abs(p) < EPS:
            if q < 0.0:
                return None
            continue
        r = q / p
        if p < 0.0:
            u0 = max(u0, r)
        else:
            u1 = min(u1, r)
        if u0 > u1:
            return None

    clipped0 = Point(p0.x + u0 * dx, p0.y + u0 * dy)
    clipped1 = Point(p0.x + u1 * dx, p0.y + u1 * dy)
    if distance(clipped0, clipped1) <= EPS:
        return None
    return clipped0, clipped1


def square_tilted_sheet(
    width: float,
    height: float,
    cell_size: float,
    theta: float,
    ligament: float,
    offset_x: float = 0.0,
    offset_y: float = 0.0,
    mirror: bool = False,
) -> Tuple[List[Polyline], List[Polyline]]:
    if width <= 0.0 or height <= 0.0:
        raise ValueError("sheet width and height must be positive")
    if cell_size <= 0.0:
        raise ValueError("cell size must be positive")

    col_min = math.floor((0.0 - offset_x) / cell_size) - 1
    col_max = math.ceil((width - offset_x) / cell_size) + 1
    row_min = math.floor((0.0 - offset_y) / cell_size) - 1
    row_max = math.ceil((height - offset_y) / cell_size) + 1

    paths: List[Polyline] = []
    for row in range(row_min, row_max + 1):
        for col in range(col_min, col_max + 1):
            origin = Point(offset_x + col * cell_size, offset_y + row * cell_size)
            for poly in square_tilted_cell_paths(origin, cell_size, theta, ligament, mirror=False):
                if len(poly.points) != 2:
                    continue
                clipped = clip_segment_to_rectangle(poly.points[0], poly.points[1], width, height)
                if clipped is not None:
                    if mirror:
                        clipped = (Point(width - clipped[0].x, clipped[0].y), Point(width - clipped[1].x, clipped[1].y))
                    paths.append(Polyline(clipped, layer=poly.layer))

    outline = [
        Polyline(
            (
                Point(0.0, 0.0),
                Point(width, 0.0),
                Point(width, height),
                Point(0.0, height),
                Point(0.0, 0.0),
            ),
            layer="OUTLINE",
        )
    ]
    return paths, outline


def triangular_tilted_sheet(
    width: float,
    height: float,
    cell_size: float,
    theta: float,
    ligament: float,
    offset_x: float = 0.0,
    offset_y: float = 0.0,
    mirror: bool = False,
    include_down_blocks: bool = False,
) -> Tuple[List[Polyline], List[Polyline], PatternParameters]:
    if width <= 0.0 or height <= 0.0:
        raise ValueError("sheet width and height must be positive")
    if cell_size <= 0.0:
        raise ValueError("cell size must be positive")

    a_over_l = math.cos(theta) - SQRT3 * math.sin(theta)
    params = PatternParameters("tilted", cell_size, a_over_l, theta=theta)
    e1 = Point(cell_size, 0.0)
    e2 = Point(0.5 * cell_size, 0.5 * SQRT3 * cell_size)

    # Conservative index bounds; generated segments are clipped to the sheet.
    col_min = math.floor((0.0 - offset_x - cell_size) / cell_size) - 2
    col_max = math.ceil((width - offset_x + cell_size) / cell_size) + 2
    row_step = 0.5 * SQRT3 * cell_size
    row_min = math.floor((0.0 - offset_y - cell_size) / row_step) - 2
    row_max = math.ceil((height - offset_y + cell_size) / row_step) + 2

    paths: List[Polyline] = []
    for row in range(row_min, row_max + 1):
        for col in range(col_min, col_max + 1):
            origin = Point(offset_x, offset_y) + e1 * col + e2 * row
            triangles = [(origin, origin + e1, origin + e2)]
            if include_down_blocks:
                triangles.append((origin + e1, origin + e1 + e2, origin + e2))
            for triangle in triangles:
                for poly in tilted_selected_paths(triangle, params, ligament):
                    if len(poly.points) != 2:
                        continue
                    clipped = clip_segment_to_rectangle(poly.points[0], poly.points[1], width, height)
                    if clipped is not None:
                        if mirror:
                            clipped = (Point(width - clipped[0].x, clipped[0].y), Point(width - clipped[1].x, clipped[1].y))
                        paths.append(Polyline(clipped, layer=poly.layer))

    outline = [
        Polyline(
            (
                Point(0.0, 0.0),
                Point(width, 0.0),
                Point(width, height),
                Point(0.0, height),
                Point(0.0, 0.0),
            ),
            layer="OUTLINE",
        )
    ]
    return paths, outline, params


def reference_template_sheet(
    width: float,
    height: float,
    cell_size: float,
    theta: float,
    ligament: float,
    offset_x: float = 0.0,
    offset_y: float = 0.0,
    mirror: bool = False,
) -> Tuple[List[Polyline], List[Polyline]]:
    """Template matching the supplied Auxetics_ULB_Generated.svg structure."""
    if width <= 0.0 or height <= 0.0:
        raise ValueError("sheet width and height must be positive")
    if cell_size <= 0.0:
        raise ValueError("cell size must be positive")

    l = cell_size
    cx = 0.5 * l
    cy = SQRT3 * l / 6.0
    x_step = 0.5 * l
    y_step = 0.5 * SQRT3 * l
    # Exact equations embedded in the ULB/FabAcademy generator.
    t_cal = ligament / l
    x = (
        (math.sin(theta - math.pi / 3.0) / math.sin(math.pi / 3.0))
        * (t_cal * math.cos(math.pi / 3.0 - theta) + (t_cal - 1.0) * math.cos(theta))
        - t_cal * math.sin(math.pi / 3.0 - theta)
    )
    y = math.tan(theta) * x + t_cal * (
        math.sin(math.pi / 3.0) - math.tan(theta) * math.cos(math.pi / 3.0)
    )
    base_start = Point(ligament * math.cos(math.pi / 3.0), ligament * math.sin(math.pi / 3.0))
    base_end = Point(x * l, y * l)

    Matrix = Tuple[Tuple[float, float, float], Tuple[float, float, float], Tuple[float, float, float]]

    def matmul(a: Matrix, b: Matrix) -> Matrix:
        return tuple(
            tuple(sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3))
            for i in range(3)
        )  # type: ignore[return-value]

    def translate(dx: float, dy: float) -> Matrix:
        return ((1.0, 0.0, dx), (0.0, 1.0, dy), (0.0, 0.0, 1.0))

    def scale(sx: float, sy: float) -> Matrix:
        return ((sx, 0.0, 0.0), (0.0, sy, 0.0), (0.0, 0.0, 1.0))

    def rotate(angle_deg: float, about_x: float = 0.0, about_y: float = 0.0) -> Matrix:
        a = math.radians(angle_deg)
        c = math.cos(a)
        s = math.sin(a)
        r = ((c, -s, 0.0), (s, c, 0.0), (0.0, 0.0, 1.0))
        return matmul(translate(about_x, about_y), matmul(r, translate(-about_x, -about_y)))

    def combine(transforms: Sequence[Matrix]) -> Matrix:
        matrix: Matrix = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
        for transform in transforms:
            matrix = matmul(matrix, transform)
        return matrix

    def transform_point(matrix: Matrix, p: Point) -> Point:
        return Point(
            matrix[0][0] * p.x + matrix[0][1] * p.y + matrix[0][2],
            matrix[1][0] * p.x + matrix[1][1] * p.y + matrix[1][2],
        )

    def to_cad(p: Point) -> Point:
        x = width - p.x if mirror else p.x
        return Point(x, height - p.y)

    cols = max(1, int(math.ceil(width / x_step)))
    rows = max(1, int(math.ceil(height / y_step)))
    outer = translate(offset_x - 0.25 * l, offset_y)

    paths: List[Polyline] = []
    for col in range(cols):
        for row in range(rows):
            cell_transforms: List[Matrix] = [outer, translate(col * x_step, row * y_step)]
            if (col + row) % 2 == 1:
                cell_transforms.extend(
                    [
                        translate(l, cy),
                        scale(-1.0, 1.0),
                        rotate(180.0, cx, cy),
                    ]
                )
            for motif_rotation in (0.0, 120.0, 240.0):
                transforms = list(cell_transforms)
                if abs(motif_rotation) > EPS:
                    transforms.append(rotate(motif_rotation, cx, cy))
                matrix = combine(transforms)
                start = transform_point(matrix, base_start)
                end = transform_point(matrix, base_end)
                clipped = clip_segment_to_rectangle(start, end, width, height)
                if clipped is not None:
                    paths.append(Polyline((to_cad(clipped[0]), to_cad(clipped[1]))))

    outline = [
        Polyline(
            (
                Point(0.0, 0.0),
                Point(width, 0.0),
                Point(width, height),
                Point(0.0, height),
                Point(0.0, 0.0),
            ),
            layer="OUTLINE",
        )
    ]
    return paths, outline


def square_reference_template_sheet(
    width: float,
    height: float,
    cell_size: float,
    theta: float,
    ligament: float,
    offset_x: float = 0.0,
    offset_y: float = 0.0,
    mirror: bool = False,
) -> Tuple[List[Polyline], List[Polyline]]:
    """Template matching the supplied square Auxetics_ULB_Generated.svg style."""
    if width <= 0.0 or height <= 0.0:
        raise ValueError("sheet width and height must be positive")
    if cell_size <= 0.0:
        raise ValueError("cell size must be positive")
    if ligament < 0.0:
        raise ValueError("t/ligament must be non-negative")

    l = cell_size
    t = ligament
    directions = (
        Point(math.cos(theta), math.sin(theta)),
        Point(-math.sin(theta), math.cos(theta)),
        Point(-math.cos(theta), -math.sin(theta)),
        Point(math.sin(theta), -math.cos(theta)),
    )
    starts = (
        Point(0.0, t),
        Point(l - t, 0.0),
        Point(l, l - t),
        Point(t, l),
    )

    t_cal = t / l
    trigo_block = t_cal * math.sin(theta) + (t_cal - 1.0) * math.cos(theta) + t_cal
    x = l * (-math.cos(theta) * trigo_block)
    y = l * (t_cal - math.sin(theta) * trigo_block)
    base_segments: List[Tuple[Point, Point]] = [
        (Point(0.0, t), Point(x, y)),
        (Point(l - t, 0.0), Point(l - y, x)),
        (Point(l, l - t), Point(l - x, l - y)),
        (Point(t, l), Point(y, l - x)),
    ]

    def to_cad(p: Point) -> Point:
        x = width - p.x if mirror else p.x
        return Point(x, height - p.y)

    cols = max(1, int(math.ceil(width / l)))
    rows = max(1, int(math.ceil(height / l)))
    paths: List[Polyline] = []

    for col in range(cols):
        for row in range(rows):
            alternate = (col + row) % 2 == 1
            for start, end in base_segments:
                transformed = []
                for point in (start, end):
                    x = l - point.x if alternate else point.x
                    y = point.y
                    transformed.append(Point(offset_x + col * l + x, offset_y + row * l + y))
                clipped = clip_segment_to_rectangle(transformed[0], transformed[1], width, height)
                if clipped is not None:
                    paths.append(Polyline((to_cad(clipped[0]), to_cad(clipped[1]))))

    outline = [
        Polyline(
            (
                Point(0.0, 0.0),
                Point(width, 0.0),
                Point(width, height),
                Point(0.0, height),
                Point(0.0, 0.0),
            ),
            layer="OUTLINE",
        )
    ]
    return paths, outline


def building_block_paths(
    triangle: Sequence[Point],
    params: PatternParameters,
    hinge: float,
    arc_segments: int,
) -> List[Polyline]:
    if hinge < 0.0:
        raise ValueError("hinge must be non-negative")
    if params.motif == "tilted":
        return tilted_paths(triangle, params, hinge)
    if params.motif == "parallel":
        return parallel_paths(triangle, params, hinge)
    if params.motif == "circular":
        return circular_paths(triangle, params, hinge, arc_segments)
    raise ValueError("Unsupported motif: {}".format(params.motif))


def triangular_grid(
    params: PatternParameters,
    rows: int,
    cols: int,
    hinge: float,
    arc_segments: int = 24,
    include_down_blocks: bool = False,
) -> Tuple[List[Polyline], List[Polyline]]:
    if rows <= 0 or cols <= 0:
        raise ValueError("rows and cols must be positive integers")

    l = params.l
    e1 = Point(l, 0.0)
    e2 = Point(0.5 * l, 0.5 * SQRT3 * l)
    paths: List[Polyline] = []

    for row in range(rows):
        for col in range(cols):
            o = e1 * col + e2 * row
            up = (o, o + e1, o + e2)
            down = (o + e1, o + e1 + e2, o + e2)
            paths.extend(building_block_paths(up, params, hinge, arc_segments))
            if include_down_blocks:
                paths.extend(building_block_paths(down, params, hinge, arc_segments))

    outline = [
        Polyline(
            (
                Point(0.0, 0.0),
                e1 * cols,
                e1 * cols + e2 * rows,
                e2 * rows,
                Point(0.0, 0.0),
            ),
            layer="OUTLINE",
        )
    ]
    return paths, outline


def bounds(polylines: Iterable[Polyline]) -> Tuple[float, float, float, float]:
    pts = [p for poly in polylines for p in poly.points]
    if not pts:
        return 0.0, 0.0, 0.0, 0.0
    return (
        min(p.x for p in pts),
        min(p.y for p in pts),
        max(p.x for p in pts),
        max(p.y for p in pts),
    )


def write_svg(
    path: Path,
    cut_paths: Sequence[Polyline],
    outline_paths: Sequence[Polyline],
    cut_width: float,
    margin: float = 5.0,
) -> None:
    all_paths = list(cut_paths) + list(outline_paths)
    xmin, ymin, xmax, ymax = bounds(all_paths)
    width = xmax - xmin + 2.0 * margin
    height = ymax - ymin + 2.0 * margin

    def transform(p: Point) -> Tuple[float, float]:
        # SVG y-axis points downward. Flip y so the geometric model stays upright.
        return p.x - xmin + margin, ymax - p.y + margin

    def svg_polyline(poly: Polyline, stroke: str, sw: float) -> str:
        pts = ["{:.6f},{:.6f}".format(*transform(p)) for p in poly.points]
        return (
            '<polyline points="{}" fill="none" stroke="{}" '
            'stroke-width="{:.6f}" stroke-linecap="round" stroke-linejoin="round" />'
        ).format(" ".join(pts), stroke, sw)

    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<svg xmlns="http://www.w3.org/2000/svg" width="{:.6f}mm" height="{:.6f}mm" '
        'viewBox="0 0 {:.6f} {:.6f}">'.format(width, height, width, height),
        '<g id="outline">',
    ]
    lines.extend(svg_polyline(poly, "#000000", max(0.05, cut_width / 2.0)) for poly in outline_paths)
    lines.append("</g>")
    lines.append('<g id="cuts">')
    lines.extend(svg_polyline(poly, "#000000", cut_width) for poly in cut_paths)
    lines.append("</g>")
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")


def write_dxf(
    path: Path,
    cut_paths: Sequence[Polyline],
    outline_paths: Sequence[Polyline],
) -> None:
    def add_polyline(lines: List[str], poly: Polyline) -> None:
        lines.extend(["0", "LWPOLYLINE", "8", poly.layer, "90", str(len(poly.points)), "70", "0"])
        for p in poly.points:
            lines.extend(["10", "{:.9f}".format(p.x), "20", "{:.9f}".format(p.y)])

    lines = ["0", "SECTION", "2", "HEADER", "0", "ENDSEC"]
    lines.extend(["0", "SECTION", "2", "TABLES"])
    lines.extend(["0", "TABLE", "2", "LAYER", "70", "2"])
    for layer, color in (("CUT", "1"), ("OUTLINE", "8")):
        lines.extend(["0", "LAYER", "2", layer, "70", "0", "62", color, "6", "CONTINUOUS"])
    lines.extend(["0", "ENDTAB", "0", "ENDSEC"])
    lines.extend(["0", "SECTION", "2", "ENTITIES"])
    for poly in outline_paths:
        add_polyline(lines, poly)
    for poly in cut_paths:
        add_polyline(lines, poly)
    lines.extend(["0", "ENDSEC", "0", "EOF"])
    path.write_text("\n".join(lines), encoding="ascii")


def write_preview_png(
    path: Path,
    cut_paths: Sequence[Polyline],
    outline_paths: Sequence[Polyline],
    scale: float = 10.0,
    margin_px: int = 32,
    line_width: int = 4,
) -> None:
    from PIL import Image, ImageDraw

    all_paths = list(cut_paths) + list(outline_paths)
    xmin, ymin, xmax, ymax = bounds(all_paths)
    width = max(1, int(math.ceil((xmax - xmin) * scale + 2.0 * margin_px)))
    height = max(1, int(math.ceil((ymax - ymin) * scale + 2.0 * margin_px)))
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)

    def transform(p: Point) -> Tuple[float, float]:
        return (p.x - xmin) * scale + margin_px, (ymax - p.y) * scale + margin_px

    for poly in outline_paths:
        draw.line([transform(p) for p in poly.points], fill=(120, 120, 120), width=2)
    for poly in cut_paths:
        draw.line([transform(p) for p in poly.points], fill=(0, 0, 0), width=line_width)

    image.save(path)


def print_summary(
    params: PatternParameters,
    rows: int,
    cols: int,
    hinge: float,
    paths: Sequence[Polyline],
    include_down_blocks: bool,
) -> None:
    print("Triangular bistable auxetic parametric model")
    print("  motif: {}".format(params.motif))
    print("  l: {:.6g} mm".format(params.l))
    print("  a/l: {:.6g}  (a = {:.6g} mm)".format(params.a_over_l, params.a))
    if params.theta is not None:
        print("  theta: {:.6g} deg".format(math.degrees(params.theta)))
    if params.radius_over_l is not None:
        print("  R/l: {:.6g}  (R = {:.6g} mm)".format(params.radius_over_l, params.radius_over_l * params.l))
    if params.width_over_l is not None:
        print("  w/l: {:.6g}  (w = {:.6g} mm)".format(params.width_over_l, params.width_over_l * params.l))
    print("  full-expansion strain estimate: {:.6g}".format(params.full_expansion_strain))
    if hinge <= EPS:
        print("  cut paths: solid continuous centerlines (hinge gap disabled)")
    else:
        print("  uncut hinge gap parameter t: {:.6g} mm".format(hinge))
    blocks_per_cell = 2 if include_down_blocks else 1
    print(
        "  grid: {} rows x {} cols rhombic cells = {} triangular blocks".format(
            rows, cols, rows * cols * blocks_per_cell
        )
    )
    print("  cut polylines: {}".format(len(paths)))


def print_square_sheet_summary(
    width: float,
    height: float,
    cell_size: float,
    theta: float,
    ligament: float,
    offset_x: float,
    offset_y: float,
    paths: Sequence[Polyline],
) -> None:
    a_over_l = square_tilted_a_over_l(theta)
    print("Square tilted bistable auxetic laser pattern")
    print("  sheet: {:.6g} x {:.6g} mm".format(width, height))
    print("  cell size l: {:.6g} mm".format(cell_size))
    print("  theta: {:.6g} deg".format(math.degrees(theta)))
    print("  derived a/l: {:.6g}  (a = {:.6g} mm)".format(a_over_l, a_over_l * cell_size))
    print("  uncut ligament t: {:.6g} mm".format(ligament))
    print("  phase offset: ({:.6g}, {:.6g}) mm".format(offset_x, offset_y))
    print("  cut polylines: {}".format(len(paths)))


def print_triangular_sheet_summary(
    width: float,
    height: float,
    cell_size: float,
    params: PatternParameters,
    ligament: float,
    offset_x: float,
    offset_y: float,
    paths: Sequence[Polyline],
    include_down_blocks: bool,
) -> None:
    print("Triangular tilted bistable auxetic laser pattern")
    print("  sheet: {:.6g} x {:.6g} mm".format(width, height))
    print("  triangular cell size l: {:.6g} mm".format(cell_size))
    print("  theta: {:.6g} deg".format(math.degrees(params.theta or 0.0)))
    print("  derived a/l: {:.6g}  (a = {:.6g} mm)".format(params.a_over_l, params.a))
    print("  uncut ligament t: {:.6g} mm".format(ligament))
    print("  phase offset: ({:.6g}, {:.6g}) mm".format(offset_x, offset_y))
    print("  include down triangles: {}".format(include_down_blocks))
    print("  cut polylines: {}".format(len(paths)))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=("reference-template", "square-reference-template", "triangular-sheet", "square-sheet", "triangular-grid"),
        default="reference-template",
        help="reference-template reproduces the triangular ULB rectangular motif; square-reference-template reproduces the square ULB rectangular motif",
    )
    parser.add_argument("--motif", choices=("tilted", "circular", "parallel"), default="tilted")
    parser.add_argument("--l", type=float, default=20.0, help="Triangular building-block side length in mm")
    parser.add_argument("--a-over-l", type=float, default=None, help="Rotating-unit side length normalized by l")
    parser.add_argument("--theta-deg", type=float, default=None, help="Tilt angle in degrees; overrides --a-over-l for tilted")
    parser.add_argument("--radius-over-l", type=float, default=None, help="Arc radius normalized by l; overrides --a-over-l for circular")
    parser.add_argument("--width-over-l", type=float, default=None, help="Parallel offset normalized by l; overrides --a-over-l for parallel")
    parser.add_argument(
        "--hinge",
        type=float,
        default=0.0,
        help="Uncut ligament gap t in mm. Use 0 for uninterrupted solid laser cut centerlines.",
    )
    parser.add_argument("--cut-width", type=float, default=0.1, help="SVG stroke width in mm; DXF remains centerline-only")
    parser.add_argument("--rows", type=int, default=2, help="Number of rhombic unit cells along the second lattice direction")
    parser.add_argument("--cols", type=int, default=2, help="Number of rhombic unit cells along the first lattice direction")
    parser.add_argument("--arc-segments", type=int, default=24, help="Segments per pi radians for circular motifs")
    parser.add_argument("--svg", type=Path, default=Path("triangular_bam.svg"))
    parser.add_argument("--dxf", type=Path, default=Path("triangular_bam.dxf"))
    parser.add_argument("--preview-png", type=Path, default=None, help="Optional PNG preview path")
    parser.add_argument("--preview-scale", type=float, default=10.0, help="Pixels per mm for --preview-png")
    parser.add_argument("--preview-line-width", type=int, default=4, help="Cut line width in pixels for --preview-png only")
    parser.add_argument("--svg-margin", type=float, default=0.0, help="SVG viewBox margin in mm")
    parser.add_argument("--no-outline", action="store_true", help="Do not include the outer sheet outline in outputs")
    parser.add_argument(
        "--include-down-blocks",
        action="store_true",
        help="Also place mirrored down-triangle blocks in each rhombic grid cell",
    )
    parser.add_argument("--sheet-width", type=float, default=100.0, help="Rectangular sheet width in mm for --mode square-sheet")
    parser.add_argument("--sheet-height", type=float, default=150.0, help="Rectangular sheet height in mm for --mode square-sheet")
    parser.add_argument("--cell-size", type=float, default=50.0, help="Square tilted unit-cell size in mm for --mode square-sheet")
    parser.add_argument("--unit-cols", type=int, default=None, help="Optional square-sheet column count; overrides --sheet-width")
    parser.add_argument("--unit-rows", type=int, default=None, help="Optional square-sheet row count; overrides --sheet-height")
    parser.add_argument("--t", type=float, default=5.0, help="Uncut ligament length in mm for --mode square-sheet")
    parser.add_argument("--angle-deg", type=float, default=None, help="Tilt angle in degrees for --mode square-sheet")
    parser.add_argument("--offset-x", type=float, default=0.0, help="Pattern phase shift in x for --mode square-sheet")
    parser.add_argument("--offset-y", type=float, default=0.0, help="Pattern phase shift in y for --mode square-sheet")
    parser.add_argument("--mirror", action="store_true", help="Mirror the square-sheet tilted motif handedness")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    params = None
    square_theta = None
    if args.mode in {"reference-template", "square-reference-template", "square-sheet", "triangular-sheet"}:
        sheet_width = args.sheet_width
        sheet_height = args.sheet_height
        if args.unit_cols is not None:
            if args.unit_cols <= 0:
                raise ValueError("--unit-cols must be positive")
            sheet_width = args.unit_cols * args.cell_size
        if args.unit_rows is not None:
            if args.unit_rows <= 0:
                raise ValueError("--unit-rows must be positive")
            sheet_height = args.unit_rows * args.cell_size
        if args.mode == "reference-template":
            theta = math.radians(args.angle_deg if args.angle_deg is not None else 15.0)
            cut_paths, outline_paths = reference_template_sheet(
                sheet_width,
                sheet_height,
                args.cell_size,
                theta,
                args.t,
                offset_x=args.offset_x,
                offset_y=args.offset_y,
                mirror=args.mirror,
            )
        elif args.mode == "square-reference-template":
            theta = math.radians(args.angle_deg if args.angle_deg is not None else 15.0)
            cut_paths, outline_paths = square_reference_template_sheet(
                sheet_width,
                sheet_height,
                args.cell_size,
                theta,
                args.t,
                offset_x=args.offset_x,
                offset_y=args.offset_y,
                mirror=args.mirror,
            )
        elif args.mode == "square-sheet":
            if args.angle_deg is not None:
                square_theta = math.radians(args.angle_deg)
            elif args.a_over_l is not None:
                square_theta = square_tilted_theta_from_a_over_l(args.a_over_l)
            else:
                square_theta = math.radians(15.0)
            square_tilted_a_over_l(square_theta)
            cut_paths, outline_paths = square_tilted_sheet(
                sheet_width,
                sheet_height,
                args.cell_size,
                square_theta,
                args.t,
                offset_x=args.offset_x,
                offset_y=args.offset_y,
                mirror=args.mirror,
            )
        else:
            if args.angle_deg is not None:
                theta = math.radians(args.angle_deg)
            elif args.a_over_l is not None:
                theta = theta_from_a_over_l(args.a_over_l)
            else:
                theta = math.radians(15.0)
            cut_paths, outline_paths, params = triangular_tilted_sheet(
                sheet_width,
                sheet_height,
                args.cell_size,
                theta,
                args.t,
                offset_x=args.offset_x,
                offset_y=args.offset_y,
                mirror=args.mirror,
                include_down_blocks=args.include_down_blocks,
            )
    else:
        params = derive_parameters(
            args.motif,
            args.l,
            args.a_over_l,
            args.theta_deg,
            args.radius_over_l,
            args.width_over_l,
        )
        cut_paths, outline_paths = triangular_grid(
            params,
            args.rows,
            args.cols,
            args.hinge,
            args.arc_segments,
            include_down_blocks=args.include_down_blocks,
        )
    if args.no_outline:
        outline_paths = []

    if args.svg:
        write_svg(args.svg, cut_paths, outline_paths, args.cut_width, margin=args.svg_margin)
    if args.dxf:
        write_dxf(args.dxf, cut_paths, outline_paths)
    if args.preview_png:
        write_preview_png(args.preview_png, cut_paths, outline_paths, scale=args.preview_scale, line_width=args.preview_line_width)

    if args.mode == "reference-template":
        theta = math.radians(args.angle_deg if args.angle_deg is not None else 15.0)
        print("Reference rectangular laser motif")
        print("  sheet: {:.6g} x {:.6g} mm".format(sheet_width, sheet_height))
        print("  cell size l: {:.6g} mm".format(args.cell_size))
        print("  theta: {:.6g} deg".format(math.degrees(theta)))
        print("  t: {:.6g} mm".format(args.t))
        print("  phase offset: ({:.6g}, {:.6g}) mm".format(args.offset_x, args.offset_y))
        print("  cut polylines: {}".format(len(cut_paths)))
    elif args.mode == "square-reference-template":
        theta = math.radians(args.angle_deg if args.angle_deg is not None else 15.0)
        print("Square reference rectangular laser motif")
        print("  sheet: {:.6g} x {:.6g} mm".format(sheet_width, sheet_height))
        print("  cell size l: {:.6g} mm".format(args.cell_size))
        print("  theta: {:.6g} deg".format(math.degrees(theta)))
        print("  t: {:.6g} mm".format(args.t))
        print("  phase offset: ({:.6g}, {:.6g}) mm".format(args.offset_x, args.offset_y))
        print("  alternating mirrored cells: true")
        print("  cut polylines: {}".format(len(cut_paths)))
    elif args.mode == "square-sheet":
        assert square_theta is not None
        print_square_sheet_summary(
            sheet_width,
            sheet_height,
            args.cell_size,
            square_theta,
            args.t,
            args.offset_x,
            args.offset_y,
            cut_paths,
        )
    elif args.mode == "triangular-sheet":
        assert params is not None
        print_triangular_sheet_summary(
            sheet_width,
            sheet_height,
            args.cell_size,
            params,
            args.t,
            args.offset_x,
            args.offset_y,
            cut_paths,
            args.include_down_blocks,
        )
    else:
        assert params is not None
        print_summary(params, args.rows, args.cols, args.hinge, cut_paths, args.include_down_blocks)
    if args.svg:
        print("  wrote SVG: {}".format(args.svg))
    if args.dxf:
        print("  wrote DXF: {}".format(args.dxf))
    if args.preview_png:
        print("  wrote preview PNG: {}".format(args.preview_png))


if __name__ == "__main__":
    main()
