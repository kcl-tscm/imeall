import transformations as quat
import numpy as np

rot_axis  = np.array([1.,1.,1.])
rot_angle = 60.*(np.pi/180.)
rot_quaternion = quat.quaternion_about_axis(rot_angle, rot_axis)
rotm           = quat.rotation_matrix(rot_angle,rot_axis)
vector    = [1,1,-2]

print 'Rotating by: ', rot_angle, 'around ', rot_axis
print 'Rotated vector ', rotm[:3,:3].dot(vector)

