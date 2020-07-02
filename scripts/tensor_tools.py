
"""
A collection of functions to rotate vectors, moment tensors, with an emphasis on Obspy data sources

A note on coordinate systems:

Obspy focal mechanisms may include any of: principal axes, nodal planes, and moment tensors. 

>principal axes are provided with length, azimuth, plunge. 
>nodal planes are provide with strike, dip, rake. 
>moment tensors are provided a set of six components, with component names that follow teh global CMT convection. 

Moment tensor components are the components of a cartesian tensor in a right handed, Local tangent plane coordinates (LTP), in which teh unit vetcors point up-south-east USE.

Whenever tensors are given as 6-compoent Voight tensors, the obspy order (from the MomentTensor docstring) is preserved:
#[tensor.m_rr, tensor.m_tt, tensor.m_pp, tensor.m_rt, tensor.m_rp, tensor.m_tp]
#[UU, SS, EE, US, UE, SE]



In order to plot, rotate and project these objects we require a consistent right-handed, local catesian coordinate system 
Even though Obspy moment tensor component given in the USE coordinates, I prefer a ENU sytesm, as this is the intuitive way I tend to see map coordinates. 

Whenever tensors are provided as full 3x3 tensors, the ENU system is used. 

What does this mean in practice:

rotation axes must be given in ENU. 
conversion of (length) aximuth, plunge, results in ENU vector components



:license: GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""
import math
import numpy as np
import obspy


def coord_transform_matrix_2d(angle_degrees):
    
    "rotates coordinate system CCW for positive angle"
    
    rad = np.deg2rad(angle_degrees)
    TM = np.matrix(([np.cos(rad), np.sin(rad)],
                   [-1*np.sin(rad), np.cos(rad)]))
    
    return TM


def lap_to_enu(length, azimuth, plunge ):
    
    """
    Convert to length azimuth,plunge to components of vector in ENU coordinate system. 
    """
    
    
    plunge_rad = np.deg2rad(plunge)
    azimuth_rad = np.deg2rad(azimuth)
    x = length * np.cos(plunge_rad) * np.sin(azimuth_rad)
    y = length * np.cos(plunge_rad) * np.cos(azimuth_rad)
    z = -length * np.sin(plunge_rad)
    return np.array([x,y,z])



def full_from_obspy_anything(obj):
    
    """
    Takes an obspy event, focal mech, moment tensor or numpy array, 
    extracts the  6-component Voight tensor
    And returns the full 3x3 Matrix. 

    
    """
    
    if "obspy.core.event" in  str(type(obj)): #assume it's an obspy event, foc mec or moment tensor
        
        voight = voigt_from_obspy_anything(obj)
        full = full_from_voight(voight)
        
    elif type(obj) is np.ndarray:
        if len(obj.ravel() == 6):             #assume its a voight matrix in obspy order
            full = full_from_voight(obj)   
            
        elif len(obj.ravel() == 9):           #assume it's already a full moment tensor
            full = obj
            

    return full
    
    
    
    #to to SEU,
    #return np.matrix(([obs_tensor.m_tt, obs_tensor.m_tp, obs_tensor.m_rt],
    #                [obs_tensor.m_tp, obs_tensor.m_pp, obs_tensor.m_rp],
    #                [obs_tensor.m_rt, obs_tensor.m_rp, obs_tensor.m_rr]))

    #to ENU
    return np.matrix(([obs_tensor.m_pp, -obs_tensor.m_tp, obs_tensor.m_rp],
                    [-obs_tensor.m_tp, obs_tensor.m_tt, -obs_tensor.m_rt],
                    [obs_tensor.m_rp, -obs_tensor.m_rt, obs_tensor.m_rr]))




def full_from_voight(voight):
    
    """
    Takes Obspy voight order moment tensor and returns full 3x3 tensor
                                              
    """
    
    #voight order
    #mt_voight = np.array([tensor.m_rr, 
    #                      tensor.m_tt, 
    #                      tensor.m_pp, 
    #                      tensor.m_rt, 
    #                      tensor.m_rp, 
    #                      tensor.m_tp])
    
    #voight to SEU
    #fullmt = np.array(([voight[1], voight[5], voight[3]],
    #                [voight[5], voight[2], voight[4]],
    #                [voight[3], voight[4], voight[0]]))
    
    #voight to ENU
    fullmt = np.array(([voight[2], -voight[5], voight[4]],
                    [-voight[5], voight[1], -voight[3]],
                    [voight[4], -voight[3], voight[0]]))
    
    
    return fullmt
    
    
    

def voight_from_full(mt):
    
    """
    Takes the full 3x3 tensor and goes back to Obpy voight order
    Useful for plotting a full MT using beach
    """
    
    #voight order
    #[tensor.m_rr, tensor.m_tt, tensor.m_pp, tensor.m_rt, tensor.m_rp, tensor.m_tp]
    
    #SEU to voight
    #obspy_voight = np.array([mt[2, 2], mt[0, 0], mt[1, 1],
    #                    mt[0, 2], mt[1, 2], mt[0, 1] ])
    
    #ENU to voight
    obspy_voight = np.array([mt[2, 2], mt[1, 1], mt[0, 0],
                        -mt[1, 2], mt[0,2], -mt[0, 1] ])
    
    return obspy_voight



def _unpack_obspy_tensorobj_to_voight(tensor):
    
    '''
    Takes a obspy MomentTensor object and returns the numoy array of components in Voight form
    '''
    
    mt_voight = np.array([tensor.m_rr, 
                          tensor.m_tt, 
                          tensor.m_pp, 
                          tensor.m_rt, 
                          tensor.m_rp, 
                          tensor.m_tp])
    
    return(mt_voight)
    


def voigt_from_obspy_anything(obj, focmec_index = 0):
    
    '''
    Takes an obspy event, focal mech, moment tensor or numpy array, 
    returns the  6-component Voight tensor
    
    '''
    
    mt = []
    
    if type(obj) is obspy.core.event.event.Event:
        fm = obj.focal_mechanisms[focmec_index]
        mt = fm.moment_tensor
        tensor = mt.tensor
    
    elif type(obj) is obspy.core.event.source.FocalMechanism:
        mt = obj.moment_tensor
        tensor = mt.tensor
    elif type(obj) is obspy.core.event.source.MomentTensor:
        mt = obj
        tensor = mt.tensor
        
    elif type(obj) is obspy.core.event.source.Tensor:
        tensor = obj
    
    mt_voight = _unpack_obspy_tensorobj_to_voight(tensor)
    
    return mt_voight





def get_rotation_matrix(axis, angle):
    """
    Returns the rotation matrix for the specified axis and angle.
    """
    axis = map(float, axis) / np.linalg.norm(axis)
    angle = np.deg2rad(angle)

    # Use c1, c2, and c3 as shortcuts for the rotation axis.
    c1 = axis[0]
    c2 = axis[1]
    c3 = axis[2]

    # Build the rotation matrix.
    rotation_matrix = np.cos(angle) * \
        np.matrix(((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))) + \
        (1 - np.cos(angle)) * np.matrix(((c1 * c1, c1 * c2, c1 * c3),
                                         (c2 * c1, c2 * c2, c2 * c3),
                                         (c3 * c1, c3 * c2, c3 * c3))) + \
        np.sin(angle) * np.matrix(((0, -c3, c2), (c3, 0, -c1), (-c2, c1, 0)))
    return rotation_matrix






def rotate_vector_tensor(tensor, angle_degrees, axis):
    
    '''
    Rotates a cartesian vector or 2nd order tensor around an axis by a given angle
    '''
    
    R  = get_rotation_matrix(axis, angle_degrees)
    
    if len(np.array(tensor).ravel()) == 3:
        #print('rotating 3x1 vector')
        rot_tensor= np.einsum('ij,j->i', R, tensor)
        #rot_tensor= np.array(R.dot(tensor)).reshape(tensor.shape)
        
    elif len(np.array(tensor).ravel()) == 9:
        rot_tensor= np.einsum('mi,nj,ij->mn', R, R, tensor)
    return rot_tensor


def mateig( mat ):
    """
    Get the eigen solution of a 3x3 moment tensor matrix 
    given as the ptb matrix and eigen values
    """
  
    # get min, int, max eigen vals, vects in that order
    (vals, vecs) = np.linalg.eigh( mat , UPLO='U')

    # sort the eigs from low to high
    ind = np.argsort(vals)
    vals = vals[ind]
    vecs = vecs[:, ind]

    # prescribe orthogonality, b = t x p
    vecs[:,1] = np.cross( vecs[:,2] , vecs[:,0] )

    return ( vecs , vals )


def sdr2mtvect( strike, dip, rake ):
    """
    Convert strike, dip, rake (radians) into mt...
    mt = [ mrr, mtt, mpp, mtp, mrp, mrt ] where r=UP, t=S, p=E  
    
    All inputs are in degrees
    
    """

    strike= np.deg2rad(strike)
    dip = np.deg2rad(dip) 
    rake =  np.deg2rad(rake)
    
    is2 = 1./math.sqrt(2.0)

    mrr =  is2*( math.sin(2*dip) * math.sin(rake) );

    mtt = -is2*( math.sin(dip)  * math.cos(rake) * math.sin(2*strike) +     
                 math.sin(2*dip) * math.sin(rake) * math.sin(strike)**2 );

    mpp =  is2*( math.sin(dip)  * math.cos(rake) * math.sin(2*strike) -
                 math.sin(2*dip) * math.sin(rake) * math.cos(strike)**2 );

    mtp = -is2*( math.sin(dip)  * math.cos(rake) * math.cos(2*strike) + 
                 0.5*math.sin(2*dip) * math.sin(rake) * math.sin(2*strike) );

    mrp =  is2*( math.cos(dip)  * math.cos(rake) * math.sin(strike)  -
                 math.cos(2*dip) * math.sin(rake) * math.cos(strike));

    mrt = -is2*( math.cos(dip)  * math.cos(rake) * math.cos(strike)  +    
                 math.cos(2*dip) * math.sin(rake) * math.sin(strike) );

    return np.array([ mrr, mtt, mpp, mrt, mrp, mtp ] )



    

def get_projected_axis(prin_axis, azimuth_degrees):
    
    """
    Given a principal axis, as a list/array in ENU compoents, return the projection of the unit vector
    On a plane orthogonal to the azimuth_degrees
    
    the projection on an orthogonal plane, 
    The two components are the aziumth-orthogonal and depthm in a local coordinate system
    
    An additional rotation is needed to correctly plot these vector components in a cartesian system
    
    """
    #fm = ev.focal_mechanisms[fmindx]
    #taxis_plunge = fm['principal_axes'][axis]['plunge']
    #taxis_az = fm['principal_axes'][axis]['azimuth']
    #taxis_length = fm['principal_axes'][axis]['length']
    
    #get the taxis in the local ENU coordinate system
    #taxis = lap_to_enu(1, taxis_az,  taxis_plunge)
    
    #Rotate the coordinate system from ENU to E'N'U such that
    #N' is parallel to the azimuth (remember azimuth is measured CW from N)
    #The E' component (TaxisNormal[0]) provides the horizontal projecion of the axis 
    #on the plane orthogonal to the azimuth (in a CW sense)
    #We then retun the projected vector E'U. \
    #Thsi could be doen better, by avoiding the 2d rotation and using the transpose of a 3d rotation matrix 
    
    
    

    TM1 = get_rotation_matrix([0,0,1], -1*azimuth_degrees).T #transform of a rot. matrix is the coord. transformation 
    rotTaxis = list(np.einsum("ij,j ->i", TM1, prin_axis))
    TaxisNormal = [rotTaxis[0], rotTaxis[-1]]
    
    
    return TaxisNormal[0], 1.*TaxisNormal[-1]


def dev_and_DC(MT):

    
    """
    Pass in any MT representation that can be handled by full_from_obspy_anything
    
    Return deviatoiric and double-couple components of the moment tensor
    
    """
    
    M = full_from_obspy_anything(MT)
    
    # isotropic part
    M_iso = np.diag(np.array([1. / 3 * np.trace(M),
                             1. / 3 * np.trace(M),
                             1. / 3 * np.trace(M)]))

    M0_iso = abs(1. / 3 * np.trace(M))

    # deviatoric part
    M_devi = M - M_iso
    
    
    #DC part
    #eigenvalues and -vectors
    eigenwtot, eigenvtot = np.linalg.eig(M_devi)

    # eigenvalues and -vectors of the deviatoric part
    eigenw1, eigenv1 = np.linalg.eig(M_devi)

    # eigenvalues in ascending order:
    eigenw = np.real(np.take(eigenw1, np.argsort(abs(eigenwtot))))
    eigenv = np.real(np.take(eigenv1, np.argsort(abs(eigenwtot)), 1))

    # eigenvalues in ascending order in absolute value!!:
    eigenw_devi = np.real(np.take(eigenw1, np.argsort(abs(eigenw1))))
    #eigenv_devi = np.real(np.take(eigenv1, np.argsort(abs(eigenw1)), 1))

    M0_devi = max(abs(eigenw_devi))

    # named according to Jost & Herrmann:
    #a1 = eigenv[:, 0]
    a2 = eigenv[:, 1]
    a3 = eigenv[:, 2]

    epsilon = 1e-13
    # if only isotropic part exists:
    if M0_devi < epsilon:
        F = 0.5
    else:
        F = -eigenw_devi[0] / eigenw_devi[2]

    M_DC = np.matrix(np.zeros((9), float)).reshape(3, 3)
    M_CLVD = np.matrix(np.zeros((9), float)).reshape(3, 3)

    M_DC = eigenw[2] * (1 - 2 * F) * (np.outer(a3, a3) - np.outer(a2, a2))
    
    return M_devi,  M_DC


def get_axis_from_obs(cat, whichaxis = 't_axis'):
    
    '''
    Pull out the t or p axis from an obspy catalog, using either the provide axis, the moment tensor,  or the strike-dip-rake
    '''
    
    axisfoundlist = []   # 1 if we could produce a principal axis, or 0 if not
    
    axislist = []        #each item is list of principal axis vector components (E,N,U) or empty of no axis could be created
    
    for i in range(len(cat)):
        
        #print(i)
        fmlist = cat[i].focal_mechanisms
        for j in range(len(fmlist)):
                
            #default if nothing is found
            axis_list_ = []
            axisfoundlist_ = 0
            
            fm = fmlist[j]
            #Check for an t-axis
            if fm['principal_axes'] is not None:
                
                taxis_plunge = fm['principal_axes'][whichaxis]['plunge']
                taxis_az = fm['principal_axes'][whichaxis]['azimuth']
                taxis_length = fm['principal_axes'][whichaxis]['length']
                
                #get the taxis in the local ENU coordinate system
                taxis = lap_to_enu(1, taxis_az,  taxis_plunge)
                    
                axislist_ = taxis
                axisfoundlist_ = 1
                break
                    
            
            #Check for a moment tensor
            if fm['moment_tensor'] is not None:
                if fm['moment_tensor']['tensor'] is not None:
                    matfull = full_from_obspy_anything(fm)
                    (V, D) = mateig(matfull)
                    paxis, baxis, taxis = V[:,0], V[:,1], V[:,2] # phi

                    if 't'in whichaxis:
                        axislist_ = taxis
                        axisfoundlist_ = 1
                    elif 'p' in whichaxis:
                        axislist_ = taxis
                        axisfoundlist_ = 1
                    else:
                        break


                    break
            
            #check for nodal planes    
            if fm['nodal_planes'] is not None:
                if fm['nodal_planes']['nodal_plane_1'] is not None:
                    no_pl = fm['nodal_planes']['nodal_plane_1']
                    voightfromsdr = sdr2mtvect(no_pl.strike, 
                                      no_pl.dip, 
                                      no_pl.rake)
                    matfull = full_from_voight(voightfromsdr)
                    (V, D) = mateig(matfull)
                    paxis, baxis, taxis = V[:,0], V[:,1], V[:,2] # phi

                    if 't'in whichaxis:
                        axislist_ = taxis
                        axisfoundlist_ = 1
                    elif 'p' in whichaxis:
                        axislist_ = taxis
                        axisfoundlist_ = 1
                    else:
                        break


                    break
                
        
        
        #Append whatever we found
        axislist.append(axislist_)       
        axisfoundlist.append(axisfoundlist_)
                
    return axislist, axisfoundlist 
    
    
def get_focmec_from_obs(cat):
    
    '''
    Pull out the moment tensor directly, or a moment tensor built from the strike-dip-rake
    '''
    
    focmeclist = []   # 1 of we could produce a principal axis, or 0 if not
    focmecfoundlist = []        #
    
    for i in range(len(cat)):
        
        #print(i)
        fmlist = cat[i].focal_mechanisms
        
        #default if nothing is found
        focmeclist_ = []
        focmecfoundlist_ = 0
        
        for j in range(len(fmlist)):
            
            fm = fmlist[j]                 
            
            #Check for a moment tensor
            if fm['moment_tensor'] is not None:
                if fm['moment_tensor']['tensor'] is not None:
                    mat = voigt_from_obspy_anything(fm)
                    focmeclist_ = mat
                    focmecfoundlist_ = 1

                    break
            
            #check for nodal planes    
            if fm['nodal_planes'] is not None:
                if fm['nodal_planes']['nodal_plane_1'] is not None:
                    no_pl = fm['nodal_planes']['nodal_plane_1']
                    mat = sdr2mtvect(no_pl.strike, 
                                      no_pl.dip, 
                                      no_pl.rake)
                    
                    

                    focmeclist_ = mat
                    focmecfoundlist_ = 1


                    break
                
        
        
        #Append whatever we found
        focmeclist.append(focmeclist_)       
        focmecfoundlist.append(focmecfoundlist_)
                
    return focmeclist, focmecfoundlist
