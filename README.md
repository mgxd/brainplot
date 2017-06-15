# Brain Plots

#### Requirements
* Mayavi
* VTK
* numpy
* matplotlib
* [Xvfbwrapper](https://github.com/cgoldberg/xvfbwrapper)
* [Nibabel 2.2.0-dev](https://github.com/nipy/nibabel/archive/master.zip)
* `git clone` this repository
* HCP 500 `rfMRI_REST1_LR_Atlas.nii`

##### Arguments [defaults in bold]
  * in_stat, path of statistic to map
  * -o: outfile, *'`in_stat`.png'*
  * -c: path to Conte69 atlas, *current working directory*
  * -r: path to HCP 500 resting atlas
  * -t: set threshold; *2.3*
  * -s: set image size; 1-smallest, 5-largest, *2*
  * -v: view of the brain: *lat*, sup, or inf
  * -i: inflation level of generated brain mesh: low, *medium*, or high
  * -w: run windowed

##### Sample Use
```python
python plotting.py /samples/faces_gt_objects.nii
```
<img src=doc/sample/lateral_quad_face_gt_object.png width=500px>

```python
python plotting.py /samples/faces_gt_objects.nii -v inf
```
<img src=doc/sample/inferior_face_gt_object.png width=500px>



###### TODO
[:)](https://github.com/mgxd/brainplot/projects/1)
