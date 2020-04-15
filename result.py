import matplotlib.pyplot as plt
import csv
import numpy as np
import math

def transform_matrix():
    yaw = np.deg2rad(124.281)
    pitch = np.deg2rad(3.4)
    roll = np.deg2rad(0.6)
    translation = np.array((-0.1, 0.25, -0.8)).reshape((3,1))

    yawMatrix = np.matrix([[math.cos(yaw), -math.sin(yaw), 0],[math.sin(yaw), math.cos(yaw), 0],[0, 0, 1]])
    pitchMatrix = np.matrix([[math.cos(pitch), 0, math.sin(pitch)],[0, 1, 0],[-math.sin(pitch), 0, math.cos(pitch)]])
    rollMatrix = np.matrix([[1, 0, 0],[0, math.cos(roll), -math.sin(roll)],[0, math.sin(roll), math.cos(roll)]])

    R = yawMatrix * pitchMatrix * rollMatrix

    T = np.concatenate((R, translation), axis=1)
    T = np.row_stack((T, [0,0,0,1]))

    return T

def transform_points(lr_target):
    T = transform_matrix()
    lr_target_transformed = lr_target.copy()
    for i in range(0, len(lr_target)):
        p = np.array([lr_target[i,1], lr_target[i,2], lr_target[i,3], 1]).reshape((4,1))
        p = np.dot(np.linalg.inv(T), p)
        lr_target_transformed[i][1:4] = p[0:3].reshape((3,))
    return lr_target_transformed

def cartesian_to_spherical(p):
    r = np.linalg.norm(p[1:4])
    azimuth = np.arctan2(p[2], p[1])
    azimuth = np.rad2deg(azimuth)
    return np.array([p[0], azimuth, r])

if __name__=='__main__':
    #read data
    file = open('data_out.csv', 'r')
    file_g = open('ground_truth.csv', 'r')
    csvCursor = csv.reader(file)
    csvCursor_g = csv.reader(file_g)
    i = 0    
    for row in csvCursor:
        if i == 0:
            row = [float(x) for x in row[0:7]]
            data = np.array(row[0:7])
        elif i > 0:
            row = [float(x) for x in row[0:7]]
            data = np.row_stack((data, row[0:7]))
        i = i + 1
    i = 0    
    for row in csvCursor_g:
        if i == 0:
            row = [float(x) for x in row[0:4]]
            ground_truth = np.array(row[0:4])
        elif i > 0:
            row = [float(x) for x in row[0:4]]
            ground_truth = np.row_stack((ground_truth, row[0:4]))
        i = i+1

    a_raw = np.rad2deg(np.arctan2(data[:,2], data[:,1]))
    a_pred = np.rad2deg(np.arctan2(data[:,4], data[:,3]))
    a_update = np.rad2deg(np.arctan2(data[:,6], data[:,5]))

    x_raw = data[:,1]
    y_raw = data[:,2]
    x_pred = data[:,3]
    y_pred = data[:,4]
    x_update = data[:,5]
    y_update = data[:,6]
    t = data[:, 0] - data[0,0]

    ground_truth = transform_points(ground_truth)

    for i in range(0, len(ground_truth)):
        if i == 0:
            s = cartesian_to_spherical(ground_truth[i][0:4])
        else:
            temp = cartesian_to_spherical(ground_truth[i][0:4])
            s = np.row_stack((s, temp))
    x_g = np.multiply(s[:,2], np.cos(s[:,1]*math.pi/180));
    y_g = np.multiply(s[:,2], np.sin(s[:,1]*math.pi/180));
    a_g = s[:,1]
    t_g = s[:,0] - data[0,0] - 0.05


    # time v.s. azimuth
    plt.figure()
    plt.scatter(t, a_raw, s = 0.2, color = 'red', label = 'raw')
    plt.scatter(t, a_pred, s = 0.2, color = 'blue', label = 'pred')
    plt.scatter(t, a_update, s = 0.2, color = 'green', label = 'update')
    plt.plot(t_g, a_g, linewidth=0.5, color = 'black', linestyle="-", label = 'ground truth')
    plt.xlabel('time')
    plt.ylabel('azimuth(deg)')
    plt.legend(loc='upper left')

    plt.savefig('ta.png', dpi = 700)
    plt.show()


    # time v.s. X
    plt.figure()
    plt.scatter(t, x_raw, s = 0.3, color = 'red', label = 'raw')
    plt.scatter(t, x_pred, s = 0.3, color = 'blue', label = 'pred')
    plt.scatter(t, x_update, s = 0.3, color = 'green', label = 'update')
    plt.plot(t_g, x_g, linewidth=0.5, color = 'black', linestyle="-", label = 'ground truth')
    plt.xlabel('time')
    plt.ylabel('x(m)')
    plt.legend(loc='upper right')

    plt.savefig('tx.png', dpi = 700)
    plt.show()


    # time v.s. y
    plt.figure()
    plt.scatter(t, y_raw, s = 0.3, color = 'red', label = 'raw')
    plt.scatter(t, y_pred, s = 0.3, color = 'blue', label = 'pred')
    plt.scatter(t, y_update, s = 0.3, color = 'green', label = 'update')
    plt.plot(t_g, y_g, linewidth=0.5, color = 'black', linestyle="-", label = 'ground truth')
    plt.xlabel('time')
    plt.ylabel('y(m)')
    plt.legend(loc='upper left')

    plt.savefig('ty.png', dpi = 700)
    plt.show()


    # time v.s. y
    plt.figure()
    plt.scatter(x_raw, y_raw, s = 0.3, color = 'red', label = 'raw')
    plt.scatter(x_pred, y_pred, s = 0.3, color = 'blue', label = 'pred')
    plt.scatter(x_update, y_update, s = 0.3, color = 'green', label = 'update')
    plt.plot(x_g, y_g, linewidth=0.5, color = 'black', linestyle="-", label = 'ground truth')
    plt.xlabel('x(m)')
    plt.ylabel('y(m)')
    plt.legend(loc='upper right')

    plt.savefig('xy.png', dpi = 700)
    plt.show()

