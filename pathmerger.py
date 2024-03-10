# Merges path elements in SVG files based on endpoint closeness considering group membership.
# Copyright (C) 2024  Stefan Schuhbäck

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import annotations
from abc import ABC, abstractmethod
from argparse import ArgumentParser
import copy
from dataclasses import dataclass, field
from typing import List, Tuple
import inkex
from enum import Enum



@dataclass 
class IndexedPoint:

    p: inkex.transforms.Vector2d 
    start: bool 
    path_index: DirectedPath

    @property
    def is_start_point(self):
        return self.start
    
    @property
    def is_end_point(self):
        return not self.is_start_point
    
    @property
    def x(self):
        return self.p.x
    
    @property
    def y(self):
        return self.p.y
    
    def distance(self, other: IndexedPoint):
        return (self.p - other.p).length()

    def is_equal_position(self, other: IndexedPoint, e:float ):
        return self.x_equal(other, e) and self.y_equal(other, e)

    
    def x_equal(self, other: IndexedPoint, e: float ):
        return abs(self.x - other.x) < e
    
    def y_equal(self, other: IndexedPoint, e: float ):
        return abs(self.y - other.y) < e
    
    def __hash__(self):
            return hash((self.x, self.y, self.start))

    def __eq__(self,other):
        return self.x == other.x and self.y == other.y and self.start == other.start
    
    @property
    def path_id(self):
        return self.path_index.path_id

    @property
    def path_id_no_dir(self):
        return self.path_index.path_id_no_dir
    
    def __repr__(self) -> str:
        if self.start:
            path = f"*{self.path_index.ascii_str()}"
        else:
            path = f"{self.path_index.ascii_str()}*"
        return f"IndexedPoint(p={self.p}, start={self.start}, path={path})"

class APath(ABC):
    """Some structure wrapping a svg path element or a list of such elements destined to be combined """

    def reverse(self) -> APath:
        """Reverse structure. Will return a copy"""
        ...

    def attach_at_start(self, other:DirectedPath|DirectedPathSegment) -> DirectedPathSegment:
        """Extends the structure and returns a copy of the combination as an DirectedPathSegment"""
        ...
    
    def attach_at_end(self, other: DirectedPath|DirectedPathSegment) -> DirectedPathSegment:
        """Extends the structure and returns a copy of the combination as an DirectedPathSegment"""
        ...

    @property
    @abstractmethod
    def start_point(self) -> IndexedPoint:
        """The first coordinate of the first element of the structure"""
        ...
    
    @property
    @abstractmethod
    def end_point(self) -> IndexedPoint:
        """The last coordinate of the last element of the structure"""
        ...


class DirectedPath(APath):
    """Wrapper class around a inkex.paths.Path element allowing reversing, extending. Each 
    extension does not modified the underling inkex.paths.Path object. Reversing will copy the object 
    and call the reverse on the inekx.paths.Path object to correctly label start and end points.
    """
    
    def __init__(self, path: inkex.paths.Path, path_id) -> None:
        self.path: inkex.paths.Path = path
        nodes = list(self.path.end_points)
        self._start: IndexedPoint = IndexedPoint(nodes[0], start=True, path_index=self)
        self._end : IndexedPoint = IndexedPoint(nodes[-1], start=False, path_index=self) 
        self.path_id: str = path_id
    
    def reverse(self) -> DirectedPath:
        """Reverse structure. Will return a copy"""
        if self.path_id.endswith("'"):
            _id = self.path_id[:-1]
        else:
            _id = self.path_id + "'"
        return DirectedPath(self.path.reverse(), path_id=_id)

    def attach_at_start(self, other:DirectedPath) -> DirectedPathSegment:
        """Extends the structure and returns a copy of the combination as an DirectedPathSegment"""
        if isinstance(other, DirectedPath):
            return DirectedPathSegment(other, self)
        elif isinstance(other, DirectedPathSegment):
            return DirectedPathSegment(*other.path_segments, self)
        else:
            raise ValueError("Wrong type")

    def attach_at_end(self, other: DirectedPath|DirectedPathSegment) -> DirectedPathSegment:
        """Extends the structure and returns a copy of the combination as an DirectedPathSegment"""
        if isinstance(other, DirectedPath):
            return DirectedPathSegment(self, other)
        elif isinstance(other, DirectedPathSegment):
            return DirectedPathSegment(self, *other.path_segments)
        else:
            raise ValueError("Wrong type")
    
    @property
    def start_point(self) -> IndexedPoint:
        """The first coordinate of the first element of the structure"""
        return self._start
    
    @property
    def end_point(self) -> IndexedPoint:
        """The last coordinate of the last element of the structure"""
        return self._end
    
    @property
    def path_id_no_dir(self):
        if self.path_id.endswith("'"):
            return self.path_id[:-1]
        else:
            return self.path_id

    def ascii_str(self):
        return f"--[{self.path_id}]-->"

    def __repr__(self) -> str:
            return f"{self.ascii_str()}"
    

class Position(Enum):
    Front = 0 
    Back = 2

class PathDir(Enum):
    Normal = 0
    Reversed = 1
    NA = 3

@dataclass
class SegmentMatch:
    """Describes where `point` is located at referenced segment
    """
    point: IndexedPoint
    at_segment_pos: Position    # Location of `point` on `segment`
    segment_index: int          # Index of `segment` in DirectedSegmentSet
    other_point: IndexedPoint   # Point that shares location of `point` but is part of another path

    def attach(self, seg: DirectedPathSegment) -> DirectedPathSegment:

        if self.at_segment_pos == Position.Front:
            # path of `other_point` is attached at the front and must therefore be an  end point
            if self.other_point.is_start_point:
                _path = self.other_point.path_index.reverse()
            else:
                _path = self.other_point.path_index
            ret = seg.attach_at_start(_path)
        else:
            # Position.Back
            # path of `other_point` is attached at the end and must therefore be a start point
            if self.other_point.is_end_point: 
                _path = self.other_point.path_index.reverse()
            else:
                _path = self.other_point.path_index
            ret = seg.attach_at_end(_path)
        
        return ret

class DirectedPathSegment(APath):

    def __init__(self, *paths: DirectedPath) -> None:
        self.path_segments: List[DirectedPath] = list(paths)
        self.path_id = "-1"
        self.parent_group = None

    def get_connected_path(self):
        connected_path: inkex.paths.Path = copy.copy(self.path_segments[0].path)
        for p in self.path_segments[1:]:
            for cmd in p.path[1:]:
                connected_path.append(cmd)

        if self.is_closed:
            connected_path.append(inkex.paths.zoneClose())
        return connected_path
            
    @property
    def start_point(self) -> IndexedPoint:
        return self.path_segments[0].start_point
    
    @property
    def end_point(self) -> IndexedPoint:
        return self.path_segments[-1].end_point
    
    def apply_dir(self, dir: PathDir) -> DirectedPathSegment:
        if dir == PathDir.Reversed:
            self.reverse()
        return self
    
    def is_closed(self, e: float ) -> bool:
        return self.start_point.is_equal_position(self.end_point, e)
    
    def can_attach_back(self, other: DirectedPath, e:float ) -> Tuple[bool, PathDir]:
        """Check if other can be appended at end  self
        
        --[self]-->{o}--[other]--> or reversed other
        --[self]-->{o}--[other*]--> or 
        """

        if self.end_point.is_equal_position(other.start_point, e):
            return True, PathDir.Normal
        elif self.end_point.is_equal_position(other.end_point, e):
            return True, PathDir.Reversed

        return False, PathDir.NA 

    def can_attach_front(self, other: DirectedPath, e: float ) -> Tuple[bool, PathDir]:
        """Check if other can be appended at start of self

        --[other]-->{o}--[self]--> or reversed
        --[other']-->{o}--[self]-->  
        """
        if self.start_point.is_equal_position(other.end_point, e):
            return True, PathDir.Normal
        elif self.start_point.is_equal_position(other.start_point, e):
            return True, PathDir.Reversed
        
        return False, PathDir.NA

    
    def reverse(self):
        _segments = []
        for seg in self.path_segments:
            _segments.append(seg.reverse())
        _segments.reverse()
        return DirectedPathSegment(*_segments)

    def attach_at_end(self, other: DirectedPath|DirectedPathSegment) -> DirectedPathSegment:
        if isinstance(other, DirectedPath):
            return DirectedPathSegment(*self.path_segments,  other)
        elif isinstance(other, DirectedPathSegment):
            return DirectedPathSegment(*self.path_segments, *other.path_segments)
        else:
            raise ValueError("Wrong type")

    def attach_at_start(self, other: DirectedPath|DirectedPathSegment) -> DirectedPathSegment:
        if isinstance(other, DirectedPath):
            return DirectedPathSegment(other, *self.path_segments)
        elif isinstance(other, DirectedPathSegment):
            return DirectedPathSegment(*other.path_segments, *self.path_segments)
        else:
            raise ValueError("Wrong type")

    @property
    def front_path(self) -> DirectedPath:
        return self.path_segments[0]
    
    @property
    def front_path_id(self) -> DirectedPath:
        return self.path_segments[0].path_id
    
    @property
    def front_path_id_no_dir(self) -> DirectedPath:
        return self.path_segments[0].path_id_no_dir
    
    @property
    def back(self) -> DirectedPath:
        return self.path_segments[-1]

    @property
    def back_path_id(self) -> str:
        return self.path_segments[-1].path_id
    
    @property
    def back_path_id_no_dir(self) -> str:
        return self.path_segments[-1].path_id_no_dir
    
    def __repr__(self) -> str:
        dots = "" if len(self.path_segments) < 3 else "..."
        return f"DirectedPathSegment(len={len(self.path_segments)}, path=[{self.front_path.ascii_str()}{dots}{self.back.ascii_str()}])"
    
@dataclass
class DirectedPathSegmentSet: 
    seg_list: List[DirectedPathSegment] = field(default_factory=list, init=False)


    def add_new_segment(self, seg:DirectedPathSegment):
        self.seg_list.append(seg)


    def create_and_append_path_to_set(self, point_a: IndexedPoint, point_b: IndexedPoint):
        """ Add new path segment from two paths, touching in [point_a, point_b]
            Can be any of the following cases. Node the `*` indicates the reversed path
            --[path_a]-->{point_a, point_b} --[path_b]-->
            --[path_a*]-->{point_a, point_b} --[path_b]-->
            --[path_a]-->{point_a, point_b} --[path_b*]-->
            --[path_a*]-->{point_a, point_b} --[path_b*]-->
        """
        if point_a.is_start_point:
            path_a: DirectedPath = point_a.path_index.reverse() # wrong way?
        else:
            path_a: DirectedPath = point_a.path_index
        if point_b.is_end_point:
            path_b: DirectedPath = point_b.path_index.reverse()
        else:
            path_b: DirectedPath = point_b.path_index
        self.add_new_segment(path_a.attach_at_end(path_b))

    def find_existing_matches(self, point_a: IndexedPoint, point_b : IndexedPoint, e:float ) -> List[SegmentMatch]:
        """ Tries to match the point pair (point_a, point_b), i.e., same location, different paths, to existing  
        `DirectedPathSegments`. Possible setups: 

        1. point pair matches once: Either point_a or point_b is attachable to *one* existing segment. 
        2. point pair matches nothing: No connection possible. Point pair creates new non-connected segment
        3. point pair matches twice to same segment: Both point_a and point_b match to the *same* 
            existing segment. One at the front, the other on the back. -> This will close the path segment! 
        4. point pair matches twice to different segments: Both point_a and point_b match to one, (not the same)
            existing segment. This will lead in the combination of two existing segments, without adding a new path!

        Returns:
            List[SegmentMatch]: 
        """
        
        ret:List[SegmentMatch] = []
        for item_idx, item in enumerate(self.seg_list):
            for i, p in enumerate([point_a, point_b]):
                if item.start_point.path_id_no_dir == p.path_id_no_dir and  item.start_point.is_equal_position(p, e):
                    _r = SegmentMatch(
                        point=p, 
                        at_segment_pos=Position.Front, 
                        segment_index=item_idx, 
                        other_point= point_b if i==0 else point_a
                    )
                    ret.append(_r)

                if item.end_point.path_id_no_dir == p.path_id_no_dir and  item.end_point.is_equal_position(p, e):
                    _r = SegmentMatch(
                        point=p, 
                        at_segment_pos=Position.Back, 
                        segment_index=item_idx, 
                        other_point= point_b if i==0 else point_a
                    )
                    ret.append(_r)
        return  ret
    
    def process_segment_matches(self, matches: List[SegmentMatch]):
        if len(matches) == 0 :
            # case 2 in `find_existing_matches` 
            return
        if len(matches) > 2 or len(matches) <1: 
            raise ValueError(f"Expected 1 or 2 elements got {len(matches)}")
        
        if len(matches) == 1:
            # case 1 in `find_existing_matches` 
            # point pair is attached to one segment only and does not attach to another one either.
            # print("MatchCase 1: append new segment ot existing one")
            m: SegmentMatch = matches[0]
            seg: DirectedPathSegment = self.seg_list.pop(m.segment_index) 
            seg = m.attach(seg)
            self.seg_list.append(seg)
        else:
            m1: SegmentMatch = matches[0] 
            m2: SegmentMatch = matches[1] 
            if m1.segment_index == m2.segment_index:
                # case 3 in `find_existing_matches`  same segment
                # todo: just close segment
                # seg: DirectedPathSegment = self.seg_list.pop(m1.segment_index)
                # print("MatchCase 3: close existing segment as point pair connects the ends.")
                pass
            else:
                # case 4 in `find_existing_matches` 
                # combine segments without adding a new path
                # print("MatchCase 4: connecting existing segments without new path element.")
                seg1: DirectedPathSegment = self.seg_list[m1.segment_index]
                seg2: DirectedPathSegment = self.seg_list[m2.segment_index]
                self.seg_list.remove(seg1)
                self.seg_list.remove(seg2)
                if m1.at_segment_pos == Position.Front:
                    # reverse seg1 to attach seg2 at back
                    seg1 = seg1.reverse()
                if m2.at_segment_pos == Position.Back:
                    # reverse seg2 to attach it at back of seg1
                    seg2 = seg2.reverse()

                _paths =  list(seg1.path_segments)
                _paths.extend(seg2.path_segments)
                seg = DirectedPathSegment(*_paths)
                self.seg_list.append(seg)

    def append_segment_or_add_segment(self,  other:DirectedPathSegment, e:float ):

        for idx, seg in enumerate(self.seg_list):

            match, dir = seg.can_attach_front(other, e)     
            if match: 
                new_seg = seg.attach_at_start(other.apply_dir(dir))
                self.seg_list[idx] = new_seg
                return True
            
            match, dir = seg.can_attach_back(other, e)
            if match:
                new_seg = seg.attach_at_end(other.apply_dir(dir), e)
                self.seg_list[idx] = new_seg
                return True
            
            return False



def print_e(e: inkex.elements.BaseElement):
    _id = e.attrib.get("id", None)
    children = len(e.getchildren())
    # print(f"{e.tag_name}[{_id}] -> {e.xml_path} #{children}")

def print_list(ee: List[inkex.elements.BaseElement]):
    for e in ee:
        print_e(e)


def path_closed(path: inkex.paths.Path):
    return isinstance(path[-1], (inkex.paths.ZoneClose, inkex.paths.zoneClose))

def get_indexed_points(path: inkex.paths.Path, i):
    """Assumes open path, no zoneClose or ZoneClose in path. Will not be checked!"""
    points = list(path.end_points)
    return [IndexedPoint(points[0], True, i), IndexedPoint(points[-1], False, i)]


def add_close_if_start_stop_same(path: inkex.paths.Path, e: float ):
    """Assumes open path, no zoneClose or ZoneClose in path. Will not be checked!"""
    if not path_closed(path):
        points: List[inkex.transforms.Vector2d] = list(path.end_points)
        diff  = points[0] - points[-1]
        if abs(diff.x) < e and abs(diff.y) < e:
            path.close()
    return path

def find_overlapping_point_pairs( points: List[IndexedPoint], e:float ) -> List[Tuple[IndexedPoint, IndexedPoint]]:

    point_pairs = []
    processed_index = set()
    i = 0
    j = i+1
    max_point = len(points)
    while j < max_point:
        if  not points[i].y_equal(points[j], e):
            # y-dist to big no other point j+1 would  overlap with i. move to next nodes
            i += 1
            j = i + 1
                
        else:
            if points[i].is_equal_position(points[j], e):
                # found match. record indices 
                if i in processed_index:
                    raise ValueError(f"Point {i} already processed. Additional match found with {j}")
                if j in processed_index:
                    raise ValueError(f"Point {j} already processed. Additional match found with {i}")
                processed_index.add(i)
                processed_index.add(j)
                point_pairs.append((points[i], points[j]))
                i += 1
                j = i + 1
            else:
                # keep incrementing j
                j += 1

    return point_pairs



class PrepareForLaserCuttingExtension(inkex.EffectExtension):

    _DEFAULT_EPSILON = 0.02

    @property
    def epsilon(self) -> float:
        return self.options.epsilon

    def add_arguments(self, pars: ArgumentParser) -> None:
        pars.add_argument("--epsilon", type=float, default= self._DEFAULT_EPSILON)
        pars.add_argument("--clearEmptyGroups", type=inkex.Boolean, default=False)

    def remove_empty_groups(self):
        empty_groups = self.svg.root.xpath("//svg:g[count(*)=0]")
        for e in empty_groups:
            if e.TAG == "g" and len(e.getchildren())==0:
                e.getparent().remove(e)
        # print_list(self.svg.root.xpath("//svg:g"))

    def get_sorted_path_groups_depth_first(self) -> List[inkex.elements.Group]:
        groups_with_paths = self.svg.root.xpath("//svg:g[count(child::svg:path)>0]")

        # print("#"*80)
        # print(f"found {len(groups_with_paths)} groups with path elements")
        
        groups_depth_first = [(p, p.xml_path.split("/")) for p in groups_with_paths]
        groups_depth_first.sort(key=lambda x: x[1], reverse=True)
        return [i[0] for i in groups_depth_first]

    @staticmethod
    def create_open_directed_path_list(path_group: inkex.elements.Group):
        # print("process ", end="")
        # print_e(path_group)
        open_paths:DirectedPath = []
        c = 0
        for child in path_group:
            if isinstance(child, inkex.elements.PathElement):
                _p = child.path.transform(child.transform)
                c += 1
                if not path_closed(_p):
                    _p = _p.to_relative()
                    open_paths.append(DirectedPath(_p, path_id=child.attrib["id"]))
                    # print("open path. use ", end="")
                    # print_e(child)
                else:
                    # print("closed path. ignore ", end="")
                    # print_e(child)
                    pass
        return open_paths, c

    @staticmethod
    def get_sorted_path_end_points(paths: List[DirectedPath]) -> List[IndexedPoint]:
        indexed_points: List[IndexedPoint] = []
        for p in paths:
            indexed_points.append(p.start_point)
            indexed_points.append(p.end_point)

        indexed_points.sort(key= lambda point: point.y, reverse=True) 
        return indexed_points

    def create_merged_segments(self, groups: List[inkex.elements.Group]):
        merged_paths: List[DirectedPathSegment] = []
        for path_group in groups:

            open_paths:List[DirectedPath]
            num_path_elements: int
            open_paths, num_path_elements = self.create_open_directed_path_list(path_group)
            # print(f"from {num_path_elements} path children found {len(open_paths)} open subpaths")

            # open_paths = [open_paths[36],open_paths[35],open_paths[34]]
            path_end_points = self.get_sorted_path_end_points(open_paths)
            
            overlapping_point_set = find_overlapping_point_pairs(path_end_points, e=self.epsilon)

            # print(f"found {len(overlapping_point_set)} overlapping points out of {len(path_end_points)} points")
            # print("#"*80)

            processed_paths: DirectedPathSegmentSet = DirectedPathSegmentSet()

            for idx, (point_a, point_b)  in enumerate(overlapping_point_set):
                # print(f"process pair {idx}: ({point_a}, {point_b})") 
                matched_segments: List[SegmentMatch]= processed_paths.find_existing_matches(point_a, point_b, e=self.epsilon)
                if len(matched_segments) > 0:
                    # append to existing segment
                    # 1. find segment to which the DirectedPath can be attached to. 
                    # 2. check if the new segment can be attached to another segment. 
                    #     If yes  do so. This happens if the new path fills a gap <----->< * new one * ><------>
                    # print(f"process pair {idx}: found matching {len(matched_segments)} segment(s) in processed_paths. Append to it/them.") 
                    processed_paths.process_segment_matches(matched_segments)

                else: 
                    # print(f"process pair {idx}: No matching segment in processed_paths found. Add new segment") 
                    processed_paths.create_and_append_path_to_set(point_a, point_b)

            for idx, seg in enumerate(processed_paths.seg_list):
                seg.path_id = f"{path_group.eid}-{idx}"
                seg.parent_group = path_group
                merged_paths.append(seg)
        return merged_paths
    
    def by_id(self, id, must_exist: bool = False) -> inkex.elements.BaseElement:
        ret = self.svg.root.xpath(f"//*[@id='{id}']")
        if len(ret) == 0:
            if must_exist:
                raise ValueError(f"found not element with id '{id}'")
            else:
                return None
        elif len(ret) > 1:
            raise ValueError(f"found {len(ret)} elements with id '{id}'. There should only be one!")
        else:
            return ret[0]
    
    def apply_merges_to_svg(self, merged_segments: List[DirectedPathSegment]):
        for seg in merged_segments:
            _parent = self.by_id(seg.parent_group.eid, must_exist=True)
            _style = None
            ns_map = None
            for path in seg.path_segments:
                _el = self.by_id(path.path_id_no_dir, must_exist=True)
                if _style is None:
                    _style = _el.attrib["style"]
                    ns_map = _el.nsmap
                _parent.remove(_el)
            
            new_elm = inkex.PathElement(attrib={"style":_style}, nsmap=ns_map)
            new_elm.set_path(seg.get_connected_path())
            new_elm.set_id(seg.path_id)
            _parent.add(new_elm)
            # delete old path segments in group

    def effect(self):
        if self.options.clearEmptyGroups:
            self.remove_empty_groups()

        groups_by_depth: List[inkex.elements.Group] = self.get_sorted_path_groups_depth_first()
        # groups_by_depth = [groups_by_depth[1]]

        merged = self.create_merged_segments(groups_by_depth)
        self.apply_merges_to_svg(merged)
        # start with inner most group 
        # get path from group, make it absolute, break ist apart to find path s

        

        
if __name__ == '__main__':
    # inkex.utils.debug(sys.argv)
    PrepareForLaserCuttingExtension().run()        