'''
Get coordinates for error ellipses, based on coordinates
x-error, y-error (1 sigma) and error correllation, and scaling 
factor (defaults to 95% confidence).
'''

import numpy as np
import math

def calc_ellipse_params(x_err, y_err, rho, scale=2.4477):
    '''
    Calculate major and minor axes and inclination of ellipse.
    returns x_size, y_size, theta.
    '''
    xy_err = rho * x_err * y_err
    covmat = np.matrix([[x_err**2,xy_err],
                        [xy_err, y_err**2]])

    eig = np.linalg.eigvals(covmat)

    theta  = 1/2 * math.atan((2*xy_err)/(x_err**2-y_err**2))
    x_size = eig[0]**0.5 * scale
    y_size = eig[1]**0.5 * scale
    if x_err >= y_err:
        theta = -theta
        
    return (x_size, y_size, theta)


def ellipse_formula(x,y,a,b,theta):
    x_t = lambda t: x + a*math.cos(t)*math.cos(theta) - b*math.sin(t)*math.sin(theta)
    y_t = lambda t: y + b*math.sin(t)*math.cos(theta) - a*math.cos(t)*math.sin(theta)
    return lambda t:[x_t(t), y_t(t)]

def ellipse(x,y,a,b,theta, num_pts=200):
#     x_t = lambda t: x + a*math.cos(t)*math.cos(theta) - b*math.sin(t)*math.sin(theta)
# #     y_t = lambda t: y + b*math.sin(t)*math.cos(theta) - a*math.cos(t)*math.sin(theta)
#     return np.array([[x_t(t) for t in np.linspace(0,2*math.pi- 2*math.pi/num_pts)],
#                      [y_t(t) for t in np.linspace(0,2*math.pi- 2*math.pi/num_pts)]]).T
    form = ellipse_formula(x,y,a,b,theta)
    return np.array([form(t) for t in np.linspace(0,2*math.pi- 2*math.pi/num_pts, num=num_pts)])

def error_ellipse_formula(x, y, x_err, y_err, rho, scale=2.4477):
    '''
    Ellipse formula generator, similar to error ellipse, but returns a function taking an angle
    (positive rotation direction) in degrees and returns the edge coordinate in tat direction.
    '''
    x_size, y_size, theta = calc_ellipse_params(x_err, y_err, rho, scale)

    return ellipse_formula(x,y,x_size,y_size,theta)

def error_ellipse(x, y, x_err, y_err, rho, scale=2.4477, num_points = 200):
    '''
    Takes data coordinates, 1 sigma errors, error correllation factor, scale
    (95% conf default (2.4477), and number of edge points).
    returns a coordinate matrix for for the error ellipse edges for each datapoint,
    in data coordinates. 
    This can be used for a matplotlib patches.
    
    Example:
    [x,y,xerr,yerr] = np.array([1,2,3]) # n-dimensional vector
    rho = np.array([0.91,0.5,.096])
    e_coords = error_ellipse(x,y,x_err/2, y_err/2, rho, num_points=20)
    e = mpl.patches.Polygon(e_coords, joinstyle='round')
    ax.add_patch(e)
    '''
    x_size, y_size, theta = calc_ellipse_params(x_err, y_err, rho, scale)
    
    return ellipse(x,y,x_size,y_size,theta,num_points)




# tests
# import matplotlib.pyplot as plt
# fig, ax = plt.subplots()
# ax.set_xlim(-8,8)
# ax.set_ylim(-8,8)
# ax.set_aspect('equal')
# e_coords = ellipse(0,0,8,8,0)
# e = mpl.patches.Polygon(e_coords)
# ax.add_patch(e)


# data = pd.DataFrame(
#        [[14.2,215],
#         [16.4,325],
#         [11.9,185],
#         [15.2,332],
#         [18.5,406],
#         [22.1,522],
#         [19.4,412],
#         [25.1,614],
#         [23.4,544],
#         [18.1,421],
#         [22.6,445],
#         [17.2,408]])



# fig,axes = plt.subplots(4,2, figsize=get_figsize('2 col',2))
# for axe, lookup in zip(axes.T,[[0,1],[1,0]]):

#     for ax,(x_mul,y_mul) in zip(axe,[[1,1],[-1,1],[1,-1],[-1,-1]]):
#         xs = data[lookup[0]]*x_mul
#         ys = data[lookup[1]]*y_mul
#         df = pd.DataFrame(np.array([xs,ys]).T)
#         rho = df.corr().loc[1,0]
#         x,y = df.mean()
#         x_err, y_err = df.std()
#         ax.plot(xs,ys,'o')
#         e_coords = error_ellipse(x,y,x_err,y_err,rho,scale=1)
#         e = mpl.patches.Polygon(e_coords,facecolor='#00000000',edgecolor='#000000FF', linewidth=1,joinstyle='round')
#         ax.annotate('{}'.format(rho),(.2,.8),xycoords='axes fraction')
#         ax.add_patch(e)


def get_bestfitdist(x, y, x_err, y_err, rho,
                    fit):
    
    shift_vec = np.matrix([[0,fit.y0]])
    l0 = np.matrix([[1,fit.alpha]])

    coords = np.array([[x,y]]) - shift_vec

    

    
    a, b, theta = calc_ellipse_params(
        x_err, y_err, rho)
    
    ellipse_form = ellipse_formula(
        x,y,a,b,theta)
    
    scale_mx_1 = np.matrix([[1,0],[0,b/a]])
    scale_mx_2 = np.matrix([[1,0],[0,a/b]])

    rot_mx = lambda t: np.matrix([[np.cos(t), - np.sin(t)],
                                  [np.sin(t), np.cos(t)]])
    
    c_pp = (coords@ rot_mx(theta) @scale_mx_1).T
    l0_pp =  (l0@ rot_mx(theta) @scale_mx_1).T
    l0_alpha = l0_pp[1,0]/l0_pp[0,0]
    n_pp = np.matrix([[-l0_pp[1,0],l0_pp[0,0]]]).T
    n_alpha = n_pp[1,0]/n_pp[0,0]
    n_0 = c_pp[1,0] - c_pp[0,0]*n_alpha
    
    x_intercept = n_0 / (l0_alpha-n_alpha)
    c_adj_pp = np.matrix([[x_intercept,l0_alpha*x_intercept]])
    
    n = ((n_pp.T * scale_mx_2) * rot_mx(-theta)) + shift_vec
    n = n/np.linalg.norm(n)
    c = ((c_pp.T * scale_mx_2) * rot_mx(-theta)) + shift_vec
    c_adj = ((c_adj_pp * scale_mx_2) * rot_mx(-theta)) + shift_vec
    theta_n = np.arccos(n[0,0])
    ell_coords = ellipse_form(theta_n)
    dir_95_conf = np.linalg.norm(np.array(c)[0]- np.array(ell_coords)) 
    misfit = np.linalg.norm(c-c_adj) * (1 if c[0,0] >= c_adj[0,0] else -1) 
    return (c_adj, misfit, n, theta_n, dir_95_conf)
    