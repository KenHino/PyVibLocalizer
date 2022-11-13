import bpy
from PIL import ImageColor
import units
import mendeleev



def plot_atom(element, x, y, z, number):
    ele = mendeleev.element(element)
    radius = ele.atomic_radius * 1.0e-02
    color = ImageColor.getcolor(ele.molcas_gv_color, "RGB")

    bpy.ops.mesh.primitive_ico_sphere_add(radius=radius, location=[x,y,z])

    name = f'{number}_{element}'
    bpy.context.object.name = name

    mat = bpy.data.materials.new(element)
    mat.diffuse_color = (*color,1)

    atom = bpy.data.object[name]
    atom.data.materials.append(mat)



if __name__ == '__main__':
    # ========= DELETE ALL mesh, light, camera, みな削除する2行 =========
    for item in bpy.data.objects:
        bpy.data.objects.remove(item)

    # ======NEW CAMERA
    bpy.ops.object.camera_add(enter_editmode=False, align='VIEW', location=(4, -6, 4), rotation=(1.3, 0, 0.5))
    ##(1.3, 0, -0.3)
    plot_atom('C', 1.0, 1.0, 1.0, 0)
    plot_atom('H', 0.0, 0.0, 0.0, 1)
