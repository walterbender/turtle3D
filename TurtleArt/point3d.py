import sys, math
import numpy as np
 
class Point3D:
    def __init__(self, x = 0, y = 0, z = 0):
        self.x, self.y, self.z = float(x), float(y), float(z)
 
    def view_calc(self, camera_position):
        #CAMERA CO-ORDINATES

        eye = np.array(camera_position)
        target = np.array([0, 0, 0])
        up = np.array([0, 1, 0])
        if eye[0] == 0 and eye[2] == 0:
            up = np.array([0, 0, -1])

        vz = eye - target
        vz = vz / np.linalg.norm(vz)
        vx = np.cross(up, vz)
        vx = vx / np.linalg.norm(vx)
        vy = np.cross(vz, vx)

        view = np.matrix([[vx[0], vx[1], vx[2], 0], \
                          [vy[0], vy[1], vy[2], 0], \
                          [vz[0], vz[1], vz[2], 0], \
                          [-np.dot(vx,eye), -np.dot(vy,eye), -np.dot(vz,eye), 1]])

        return view

    def perspective2(self, win_width, win_height):

        near = 0.001
        far = 100.0
        angle_of_view = 45.0
        aspect_ratio = win_width / win_height
   
        size = near * math.tan(math.radians(angle_of_view) / 2.0)
        left = -size
        right = size
        bottom = -size / aspect_ratio
        top = size / aspect_ratio

        proj = np.matrix([[2*near/(right - left), 0.0, (right+left)/(right-left), 0.0], \
                          [0.0, 2*near/(top - bottom), (top+bottom)/(top-bottom), 0.0], \
                          [0.0, 0.0, -1.0*(far+near)/(far-near), (-2.0*far*near)/(far-near)], \
                          [0.0, 0.0, -1.0, 0.0]])
        return proj

    def perspective(self, x, y, z, win_width, win_height):

        fov = 512
        viewer_distance = 512
        factor = fov / (viewer_distance + z)
        x = -x * factor + win_width / 2
        y = -y * factor + win_height / 2
        return [x, y, 1]

    
    def project(self, win_width, win_height, camera_position):
        """ Transforms this 3D point to 2D using a perspective projection. """

        # START CALCULATE SCREEN COORDINATES
        point = np.matrix([[self.x], \
                           [self.y], \
                           [self.z], \
                           [1.]])

        view = self.view_calc(camera_position)

        # Scale Part of Model Matrix
        scale = 1
        scalemat = np.matrix([[scale, 0, 0, 0], \
                             [0, scale, 0, 0], \
                             [0, 0, scale, 0], \
                             [0, 0, 0, 1]])

        self.modelView = view * scalemat
        self.eyeView = self.modelView * point
        final = self.perspective(self.eyeView[0,0], self.eyeView[1,0], self.eyeView[2,0], win_width, win_height)
        scrn_x = final[0]
        scrn_y = final[1]
    
        x, y = scrn_x, scrn_y
        return Point3D(x, y, 1)
