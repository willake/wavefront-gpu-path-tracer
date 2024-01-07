# Hui GPU Path Tracer

This renderer is based on a [cpu ray tracer](https://github.com/willake/cpu-ray-tracer) I made, which is extended from [template project](https://github.com/jbikker/tmpl8rt_UU) of [Dr. Jacco Bikker](https://github.com/jbikker). The following feature is developed:
**Graphics related**

-   CPU Path Tracing
-   GPU Path Tracing
-   BVH(acceleration structure)
    -   BLAS/TLAS
-   skybox
-   multiple light sources
-   reflection / refraction
-   NEE (Next Event Estimation)
-   Anti-aliasing
    **Engine related**
-   .obj model loading
-   material (reflection, refraction, texture)
-   texture loading & mapping
-   scene file loading

## IDE

Visual Studio 2022

## Path Tracer

<img src="./assets/readme/view-2-path.JPG" width="600" height="400" alt="path tracer"/>
<img src="./assets/readme/view-1-path.JPG" width="600" height="400" alt="path tracer"/>
<img src="./assets/readme/path-tracer-2.JPG" width="600" height="400" alt="path tracer"/>

## How to run

Open `tmpl_2022-rt.sln` with Visual Studio. Select a project (either CPUPathTracer or GPUPathTracer). Build.

### Inspect traversal

Toggle `Inspect Traversal` checkbox in the panel will turn into traversal debug mode.
<img src="./assets/readme/traversal.jpg" width="600" height="400" alt="path tracer"/>
<img src="./assets/readme/view-2-bvhsah-traversal.JPG" width="600" height="400" alt="path tracer"/>

## How to configure

### Scene

There are several scenes available in `assets` folder. In `renderer.h`, the user can set the path to the scene file and start the program. The scene will be loaded automatically.
A scene template looks like the following

```
<?xml version="1.0" encoding="UTF-8"?>
<scene>
    <!-- Scene Information -->
    <scene_name>wok scene</scene_name>
	<light_position>
		<x>0.0</x>
		<y>4.5</y>
		<z>2.0</z>
	</light_position>
	<plane_texture_location>../assets/textures/Stylized_Wood_basecolor.tga</plane_texture_location>
    <skydome_location>../assets/industrial_sunset_puresky_4k.hdr</skydome_location>
	<!-- Light List -->
    <lights>
        <!-- Light 0 -->
        <light>
			<position>
                <x>1.3</x>
                <y>4.5</y>
                <z>2.0</z>
            </position>
			<rotation>
                <x>0.0</x>
                <y>0.0</y>
                <z>0.0</z>
            </rotation>
			<size>1</size>
			<color>
				<x>240</x>
                <y>0</y>
                <z>0</z>
			</color>
        </light>
		<light>
			<position>
                <x>-3.3</x>
                <y>0</y>
                <z>2.8</z>
            </position>
			<rotation>
                <x>0.0</x>
                <y>0.0</y>
                <z>90.0</z>
            </rotation>
			<size>1</size>
			<color>
				<x>0</x>
                <y>0</y>
                <z>220</z>
			</color>
        </light>
        <!-- Add more lights as needed -->
    </lights>
    <!-- Object List -->
    <objects>
        <!-- Object 1 -->
        <object>
            <mesh_idx>0</mesh_idx>
            <material_idx>0</material_idx>
            <position>
                <x>0.0</x>
                <y>1</y>
                <z>2.8</z>
            </position>
            <rotation>
                <x>-90.0</x>
                <y>0.0</y>
                <z>180.0</z>
            </rotation>
			<scale>
                <x>0.005</x>
                <y>0.005</y>
                <z>0.005</z>
            </scale>
        </object>

        <!-- Add more objects as needed -->
    </objects>
	<!-- Mesh List -->
    <meshes>
        <!-- Mesh 0 -->
        <mesh>
			<id>0</id>
            <model_location>../assets/urna.obj</model_location>
			<scale>
                <x>0.005</x>
                <y>0.005</y>
                <z>0.005</z>
            </scale>
        </mesh>

        <!-- Add more objects as needed -->
    </meshes>
	<materials>
		<!-- Material 0 -->
		<material>
			<id>0</id>
			<reflectivity>0.0</reflectivity>
			<refractivity>0.0</refractivity>
			<absorption>
                <x>0.0</x>
                <y>0.0</y>
                <z>0.0</z>
            </absorption>
			<texture_location>../assets/textures/urna.jpg</texture_location>
		</material>
	</materials>
</scene>
```

## Original README.md from template

This template is intended for students of Utrecht University.

**Please refer to "\_ getting started.pdf" for instructions.**

Code is fully public domain. Use as you please.

Contact me at bikker.j@gmail.com.
