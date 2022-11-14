try:
    import bpy
except ImportError:
    print('You cannot use blender by import error')
import math
import mendeleev
import numpy as np
from PIL import ImageColor
from typing import List, Optional, Tuple, Union

import units

class Visualizer():
    def __init__(self, vibration, arrow_scale: float):
        self.freq = vibration.freq
        self.disp = vibration.disp
        self.geom = vibration.geom
        self.bond = vibration.bond
        self.coord = vibration.coord
        self.atom = vibration.atom
        self.atom_number = vibration.atom_number
        self.natom = len(self.geom)
        self.number = 0
        self.bond_number = 0

        for item in bpy.data.objects:
            bpy.data.objects.remove(item)
           
        for item in bpy.data.collections:
            if item.name != 'Collection':
                bpy.data.collections.remove(item)

        bpy.ops.object.camera_add(enter_editmode=False, align='VIEW', 
                location=(4, -6, 4), rotation=(1.3, 0, 0.5))
        self.scale = arrow_scale

    def angles(self, vec1: np.ndarray, vec2: np.ndarray
            ) -> Tuple[float, float, float]:
        v = vec2 - vec1
        r = np.linalg.norm(v)
        if abs(r) > 1.e-16:
            v /= r
        theta = math.acos(v[2])
        if theta < 0:
            theta += 2*math.pi
        phi = math.atan2(v[1], v[0])
        if phi < 0:
            phi += 2*math.pi
        return (r, theta, phi)

    def plot_atom(self, element: str, 
            xyz: Tuple[float, float, float], number: int):
        ele = mendeleev.element(element)
        x, y, z = (np.array(xyz) / units.ANGSTROM).tolist()
        radius = ele.atomic_radius * 1.0e-02
        color = ImageColor.getcolor(ele.molcas_gv_color, "RGB")

        bpy.ops.mesh.primitive_uv_sphere_add(radius=radius, location=[x, y, z],
                segments=128, ring_count=64)

        name = f'{number}_{element}'
        bpy.context.object.name = name

        mat = bpy.data.materials.new(element)
        mat.diffuse_color = (color[0]/255, color[1]/255, color[2]/255, 0.95)

        atom = bpy.data.objects[name]
        atom.data.materials.append(mat)
        self.add_collection('molecule')

        if False:#self.atom_number:
            bpy.ops.object.text_add(location=[x,y,z])
            ob = bpy.context.object
            ob.data.body = f"{number}"
            bpy.context.object.name = f"_{number}_{element}"
            mat = bpy.data.materials.new('black')
            mat.diffuse_color = (0, 0, 0, 1)

            atom = bpy.data.objects[f"_{number}_{element}"]
            atom.data.materials.append(mat)
            self.add_collection('molecule')

    def plot_bonds(self, 
        bonds: List[List[
            Union[Tuple[float, float, float], Tuple[float, float, float]]]]):
        mat = bpy.data.materials.new('bond')
        mat.diffuse_color = (1, 1, 1, 1)
        for bond in bonds:
            vec1 = np.array(bond[0]) / units.ANGSTROM
            vec2 = np.array(bond[1]) / units.ANGSTROM
            r, theta, phi = self.angles(vec1, vec2)
            bpy.ops.mesh.primitive_cylinder_add(location=tuple((vec1+vec2)/2), 
                    radius=0.05, vertices=128, depth=r, rotation=(0, theta, phi))

            self.bond_number += 1
            name = f'bond_{self.bond_number}'
            bpy.context.object.name = name
            self.add_collection('bonds')
        self.join_object('bonds')
        bpy.context.object.name = '_bonds'
        _bonds = bpy.data.objects['_bonds']
        _bonds.data.materials.append(mat)

    def plot_arrow(self, starts: List[float], vectors: List[float], 
            name: str, scale: Optional[float] =1.0):
        starts = np.array(starts) / units.ANGSTROM
        vectors = np.array(vectors).reshape(self.natom, 3) / units.ANGSTROM
        mat = bpy.data.materials.new('vec')
        mat.diffuse_color = (0, 1, 0, 1)
        for start, vector in zip(starts, vectors):
            vec1 = np.array(start)
            vec2 = vec1 + np.array(vector)*scale
            r, theta, phi = self.angles(vec1, vec2)
            self.number += 1
            if abs(r) < 1.e-16:
                continue

            bpy.ops.mesh.primitive_cylinder_add(location=tuple((vec1+vec2)/2), 
                    radius=0.10, vertices=128, depth=r, rotation=(0, theta, phi))
            name_bar = f'{name}_{self.number}_bar'
            bpy.context.object.name = name_bar

            self.add_collection(name)
            bpy.ops.mesh.primitive_cone_add(location=tuple(vec2), 
                radius1=0.25, vertices=128, depth=0.3, rotation=(0, theta, phi))

            name_cone = f'{name}_{self.number}_cone'
            bpy.context.object.name = name_cone

            self.add_collection(name)
        self.join_object(name)
        bpy.context.object.name = f'mode_{name}'
        obj = bpy.data.objects[f'mode_{name}']
        obj.data.materials.append(mat)

        self.exclude_collection(name)

    def add_collection(self, name: str):
        try:
            collection = bpy.data.collections[name]
        except KeyError:
            collection = bpy.data.collections.new(name=name)
            scene_collection = bpy.context.scene.collection
            scene_collection.children.link(collection)
        obj = bpy.context.active_object
        collection.objects.link(obj)
        collection = bpy.data.collections['Collection'].objects.unlink(obj)

    def exclude_collection(self, name: str):
        bpy.context.view_layer.layer_collection.children[name].exclude = True

    def join_object(self, name: str):
        bpy.ops.object.select_all(action='DESELECT')
        coll = bpy.data.collections[name]
        for obj in coll.objects:
            if obj.type == 'MESH':
               obj.select_set(True)
        bpy.ops.object.join()
        bpy.ops.object.select_all(action='DESELECT')

    def show(self):
        for i, atom in enumerate(self.geom):
            self.plot_atom(element=atom[0], xyz=atom[1], number=i)
        self.plot_bonds(self.bond)

        for i, f in enumerate(self.freq):
            if math.isnan(f):
                f = 0.0
            name = f'{i} : {int(f / units.CM1)} cm-1'
            self.plot_arrow(self.coord, self.disp[i], name=name, scale=self.scale)
