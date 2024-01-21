# Hui Wavefront GPU Path Tracer

<img src="./assets/readme/base-scene.JPG" width="600" height="400" alt="Overview"/>

This renderer is based on a [cpu ray tracer](https://github.com/willake/cpu-ray-tracer) I made, which is extended from [template project](https://github.com/jbikker/tmpl8rt_UU) of [Dr. Jacco Bikker](https://github.com/jbikker). The following feature is developed:

**Graphics related**
-   CPU Path Tracing
-   wavefront GPU Path Tracing
-   BVH(acceleration structure)
    -   BLAS/TLAS
-   skybox
-   multiple light sources
-   reflection / refraction
-   NEE (Next Event Estimation)
-   anti-aliasing
-   microfacet BRDF
-   gamma correction
-   
**Engine related**
-   .obj model loading
-   material (reflection, refraction, texture)
-   texture loading & mapping
-   scene file loading

## IDE

Visual Studio 2022

## Screenshots

![Hollow Knight](./assets/readme/hollow-knight.JPG)
![Bright](./assets/readme/gamma-correction.JPG)
![Depth of Field](./assets/readme/depth-of-field.JPG)
## How to run

Open `tmpl_2022-rt.sln` with Visual Studio. Select a project (either CPUPathTracer or GPUPathTracer). Build.

### Inspect traversal

Toggle `Inspect Traversal` checkbox in the panel will turn into traversal debug mode.
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
	<plane_texture_location>../assets/textures/Stylized_Wood_basecolor.tga</plane_texture_location>
    <skydome_location>../assets/harties_4k.hdr</skydome_location>
	<!-- Light List -->
    <lights>
        <!-- Light 0 -->
        <light>
	    <position>
                <x>0</x>
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
		   <x>24</x>
                   <y>24</y>
                   <z>22</z>
                </color>
        </light>
		<light>
                <position>
                   <x>3.3</x>
                   <y>0</y>
                   <z>2</z>
               </position>
               <rotation>
                  <x>0.0</x>
                  <y>120.0</y>
                  <z>90.0</z>
              </rotation>
                  <size>1</size>
                  <color>
                        <x>48</x>
                        <y>48</y>
                        <z>44</z>
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
        </object>
        <object>
			<mesh_idx>2</mesh_idx>
            <material_idx>1</material_idx>
            <position>
                <x>1.5</x>
                <y>2</y>
                <z>1</z>
            </position>
            <rotation>
                <x>0</x>
                <y>0.0</y>
                <z>0.0</z>
            </rotation>
        </object>
        <!-- Add more objects as needed -->
    </objects>
	<!-- Mesh List -->
    <meshes>
        <!-- Mesh 0 -->
        <mesh>
            <model_location>../assets/urna.obj</model_location>
			<scale>
                <x>0.005</x>
                <y>0.005</y>
                <z>0.005</z>
            </scale>
        </mesh>
		<mesh>
            <model_location>../assets/sphere.obj</model_location>
			<scale>
                <x>0.5</x>
                <y>0.5</y>
                <z>0.5</z>
            </scale>
        </mesh>
        <mesh>
            <model_location>../assets/cup.obj</model_location>
			<scale>
                <x>0.2</x>
                <y>0.2</y>
                <z>0.2</z>
            </scale>
        </mesh>
        <!-- Add more objects as needed -->
    </meshes>
	<materials>
		<!-- Material 0 -->
		<material>
			<reflectivity>0.0</reflectivity>
			<refractivity>0.0</refractivity>
			<absorption>
                <x>0.0</x>
                <y>0.0</y>
                <z>0.0</z>
            </absorption>
			<roughness>0.3</roughness>
			<texture_location>../assets/textures/urna.jpg</texture_location>
		</material>
		<material>
			<reflectivity>0.0</reflectivity>
			<refractivity>1.0</refractivity>
			<absorption>
                <x>0.3</x>
                <y>0.3</y>
                <z>0.3</z>
            </absorption>
		<texture_location>../assets/textures/white.png</texture_location>
		</material>
	</materials>
</scene>
```

## Original README.md from template

This template is intended for students of Utrecht University.

**Please refer to "\_ getting started.pdf" for instructions.**

Code is fully public domain. Use as you please.

Contact me at bikker.j@gmail.com.
