import numpy as np

def euler_to_rot(roll, pitch, yaw):
    """
    Build a rotation matrix from roll, pitch, yaw (intrinsic rotations about x, y, z).
    Angles in radians.

    Combined rotation = Rz(yaw) * Ry(pitch) * Rx(roll).
    """
    # Rotation about X-axis (roll)
    Rx = np.array([
        [1.0,           0.0,            0.0],
        [0.0,  np.cos(roll),  -np.sin(roll)],
        [0.0,  np.sin(roll),   np.cos(roll)]
    ])

    # Rotation about Y-axis (pitch)
    Ry = np.array([
        [ np.cos(pitch), 0.0, np.sin(pitch)],
        [ 0.0,           1.0,           0.0],
        [-np.sin(pitch), 0.0, np.cos(pitch)]
    ])

    # Rotation about Z-axis (yaw)
    Rz = np.array([
        [ np.cos(yaw), -np.sin(yaw), 0.0],
        [ np.sin(yaw),  np.cos(yaw), 0.0],
        [ 0.0,          0.0,         1.0]
    ])

    # Combined rotation
    return Rz @ Ry @ Rx

def ConstructProcNoise(sigma, dt):
    Q = np.zeros((6, 6))
    Q[:3, :3] = np.eye(3) * sigma**2 * dt**3 / 3
    Q[:3, 3:] = np.eye(3) * sigma**2 * dt**2 / 2
    Q[3:, :3] = np.eye(3) * sigma**2 * dt**2 / 2
    Q[3:, 3:] = np.eye(3) * sigma**2 * dt

    return Q


def rover_clock_imu_dynamics(t, current_state, imu_data, dt, imu_model, clock_model):
    """
    Dynamics function for the rover clock model with IMU data.
    """
    pos_idx = slice(0, 3)
    vel_idx = slice(3, 6)
    clk_idx = slice(6, 8)
    att_idx = slice(8, 11)
    omega_idx = slice(11, 14)
    bias_acc_imu_idx = slice(14, 17)
    bias_gyro_imu_idx = slice(17, 20)

    posvel_idx = slice(0, 6)
    attomega_idx = slice(8, 14)

    p = current_state[pos_idx]
    v = current_state[vel_idx]
    clk = current_state[clk_idx]
    att = current_state[att_idx]
    omega = current_state[omega_idx]
    bias_acc_imu = current_state[bias_acc_imu_idx]
    bias_gyro_imu = current_state[bias_gyro_imu_idx]

    # compute next clock state
    clk_next = np.zeros(2)
    clk_next[0] = clk[0] + clk[1] * dt
    clk_next[1] = clk[1]

    # estimated imu measurements
    a_imu_meas = imu_data[0:3]
    omega_imu_meas = imu_data[3:6]

    a_imu_est = a_imu_meas - bias_acc_imu
    omega_imu_est = omega_imu_meas - bias_gyro_imu

    # compute the body-to-inertial rotation matrix
    roll, pitch, yaw = att
    R = euler_to_rot(roll, pitch, yaw)

    # rotate measured accelarion
    a_inertial = R @ a_imu_est

    # compute next states
    p_next = p + v * dt + 0.5 * a_inertial * dt**2
    v_next = v + a_inertial * dt

    # compute angular velocity
    omega_next = omega + omega_imu_est * dt
    att_next = att + omega * dt + 0.5 * omega_imu_est * dt**2
    # wrap to [-pi, pi]
    for i in range(3):
        if att_next[i] > np.pi:
            att_next[i] -= 2*np.pi
        if att_next[i] < -np.pi:
            att_next[i] += 2*np.pi

    # compute next biases
    bias_acc_imu_next = bias_acc_imu
    bias_gyro_imu_next = bias_gyro_imu

    # construct next state
    next_state = np.hstack((p_next, v_next, clk_next, att_next, omega_next, bias_acc_imu_next, bias_gyro_imu_next))

    # construct stm
    stm = np.eye(20)
    stm[pos_idx, vel_idx] = np.eye(3) * dt                # position, velocity
    stm[clk_idx, clk_idx] = np.array([[1, dt], [0, 1]])   # clock
    stm[vel_idx, bias_acc_imu_idx] = -R * dt              # velocity, bias_acc_imu
    stm[att_idx, omega_idx] = np.eye(3) * dt              # attitude, omega
    stm[omega_idx, bias_gyro_imu_idx] = np.eye(3) * dt    # omega, bias_gyro_imu

    # process noise matrix
    Q = np.zeros((20, 20))
    Q[bias_acc_imu_idx, bias_acc_imu_idx] =  (imu_model.sigma_ba * np.sqrt(dt))**2 * np.eye(3)
    Q[bias_gyro_imu_idx, bias_gyro_imu_idx] = (imu_model.sigma_bg *  np.sqrt(dt))**2 * np.eye(3)
    
    sigma_a = imu_model.sigma_a / np.sqrt(dt)
    sigma_g = imu_model.sigma_g / np.sqrt(dt)

    Q[posvel_idx, posvel_idx] = ConstructProcNoise(sigma_a, dt)
    Q[attomega_idx, attomega_idx] = ConstructProcNoise(sigma_g, dt)

    Q[clk_idx, clk_idx] = clock_proc_noise(clock_model.sigma1, clock_model.sigma2, dt)

    return next_state, stm, Q


def rotation_to_rpy(R):
    """
    Extract roll-pitch-yaw (radians) from a 3x3 rotation matrix R
    assuming the sequence: Rz(yaw) * Ry(pitch) * Rx(roll).
    
    In that convention:
        yaw   = atan2(R[1,0], R[0,0])
        pitch = -asin(R[2,0])
        roll  = atan2(R[2,1], R[2,2])
    """
    # Guard against numerical issues (clamp the value in [-1, 1] for arcsin)
    # pitch = -asin(R[2,0]) is a common approach for that sequence.
    r20 = np.clip(-R[2, 0], -1.0, 1.0)
    pitch = np.arcsin(r20)

    # Yaw
    yaw = np.arctan2(R[1, 0], R[0, 0])

    # Roll
    roll = np.arctan2(R[2, 1], R[2, 2])

    return roll, pitch, yaw


class IMUModel: 
    def __init__(self, sigma_a, sigma_g, sigma_ba, sigma_bg):
        self.sigma_a = sigma_a
        self.sigma_g = sigma_g
        self.sigma_ba = sigma_ba
        self.sigma_bg = sigma_bg

class LN200S(IMUModel):
    # from the spec sheet
    # https://cdn.northropgrumman.com/-/media/Project/Northrop-Grumman/ngc/what-we-do/space/ln-200s-ln-200hps-imu/LN-200S-LN-200HPS-inertial-measurement-unit-imu-datasheet.pdf?rev=d798f713276f4b83a2f099179646d9d3&_gl=1*1nvejwu*_ga*MTQ3ODExODcxMS4xNzM3MTU5NjE4*_ga_7YV3CDX0R2*MTczNzE1OTYxOC4xLjAuMTczNzE1OTYxOC42MC4wLjA.
    # conversions
    # https://stechschulte.net/2023/10/11/imu-specs.html#:~:text=IMU%20noise%20parameter%20units&text=The%20noise%20model%20for%20both,parameterized%20by%20a%20standard%20deviation.
    def __init__(self):
        MUG = 9.80e-6
        sigma_a = 35 * MUG    # mu g /sqrt(Hz) --> m/s^2 1/sqrt(Hz)
        sigma_ba = 300 * MUG  # mu g sqrt(Hz)  --> m/s^2 sqrt(Hz)

        DEGHRSQRT = np.pi/180/60
        DEGHR = np.pi/180/3600

        sigma_g = 0.07 * DEGHRSQRT # rad/s 1/sqrt(Hz)
        sigma_bg = 1.0 * DEGHR       # rad/s sqrt(Hz)

        print(f"LN200S IMU model:")
        print(f"  Gyro White Noise : {sigma_g:.2e} rad/s 1/sqrt(Hz)")
        print(f"  Gyro Bias        : {sigma_bg:.2e} rad/s sqrt(Hz)")
        print(f"  Acc White Noise  : {sigma_a:.2e} m/s^2 1/sqrt(Hz)")
        print(f"  Acc Bias         : {sigma_ba:.2e} m/s^2 sqrt(Hz)")

        super().__init__(sigma_a, sigma_g, sigma_ba, sigma_bg)


class ClockModel:
    def __init__(self, sigma1, sigma2):
        self.sigma1 = sigma1
        self.sigma2 = sigma2


class CSAC(ClockModel):
    def __init__(self):
        sigma1 = 3e-10
        sigma2 = 3e-12
        super().__init__(sigma1, sigma2)

class MiniRafs(ClockModel):
    def __init__(self):
        sigma1 = 1e-11
        sigma2 = 1e-15
        super().__init__(sigma1, sigma2)


def clock_proc_noise(sigma1, sigma2, dt):
    Q = np.zeros((2, 2))
    Q[0, 0] = sigma1**2 * dt + sigma2**2 * dt**3 / 3
    Q[0, 1] = sigma2**2 * dt**2 / 2
    Q[1, 0] = sigma2**2 * dt**2 / 2
    Q[1, 1] = sigma2**2 * dt

    return Q


def compute_rover_state_kinematics(times, positions, normals, imu_model, clock_model):
    """
    Given:
      - positions: List/array of shape (N, 3) with rover 3D positions in an inertial frame.
      - times:     List/array of length N with time stamps (seconds).
      - normals:   List/array of shape (N, 3) with local surface normals in the same inertial frame.

    Returns:
        - state:     Array of shape (N, 18) with rover states at each time step.
        - imu_measurement: Array of shape (N, 6) with IMU measurements at each time step.
        
    The rover frame is defined so that:
      - x-axis = direction of travel (velocity direction),
      - z-axis = local surface normal (pointing 'up'),
      - y-axis = z x x, giving a right-handed frame.
    """
    pos_idx = slice(0, 3)
    vel_idx = slice(3, 6)
    clk_idx = slice(6, 8)
    att_idx = slice(8, 11)
    omega_idx = slice(11, 14)
    bias_acc_imu_idx = slice(14, 17)
    bias_gyro_imu_idx = slice(17, 20)

    # compute the velocity and acceleration at each point
    n = len(times)
    v = np.zeros((n, 3))
    a = np.zeros((n, 3))

    for i in range(1, n):
        dt = times[i] - times[i-1]
        v[i] = (positions[i] - positions[i-1]) / dt

    v[0] = v[1] - (v[2] - v[1])
    
    for i in range(1, n):
        a[i] = (v[i] - v[i-1]) / dt

    # for the first timestep assume its same as the second
    a[0] = a[1] - (a[2] - a[1])

    # compute clock
    clk_bias = 0.0
    clk_drift = 0.0
    clock_state = np.zeros((n, 2))
    clock_state[0] = np.array([clk_bias, clk_drift])
    for i in range(1, n):
        dt = times[i] - times[i-1]
        Q_clk = clock_proc_noise(clock_model.sigma1, clock_model.sigma2, dt)
        noise = np.random.multivariate_normal([0, 0], Q_clk)
        clk_bias += clk_drift * dt + noise[0]
        clk_drift += noise[1]
        clock_state[i] = np.array([clk_bias, clk_drift])

    # -- 3) For each epoch, build the local coordinate axes: x, y, z ----------
    #    x-axis = velocity direction
    #    z-axis = surface normal
    #    y-axis = z x x  (right-hand rule)
    att = np.zeros((n, 3))
    for i in range(n):
        # Current velocity
        speed = np.linalg.norm(v[i])
        x_axis = v[i] / speed

        # Surface normal
        z_axis = normals[i].copy()
        norm_z = np.linalg.norm(z_axis)
        z_axis /= norm_z

        # Construct y = z x x (to ensure right-handed orientation)
        y_axis = np.cross(z_axis, x_axis)
        norm_y = np.linalg.norm(y_axis)
        y_axis /= norm_y

        # Re-orthonormalize slightly:
        # Recompute z = y x x to ensure perfect orthonormality
        z_axis = np.cross(x_axis, y_axis)

        R = np.stack((x_axis, y_axis, z_axis), axis=1)
        roll, pitch, yaw = rotation_to_rpy(R)
        att[i, 0] = roll
        att[i, 1] = pitch
        att[i, 2] = yaw

    # compute omega (angular velocity)
    omega = np.zeros((n, 3))
    domega = np.zeros((n, 3))
    for i in range(1, n):
        dt = times[i] - times[i-1]
        att_diff = att[i] - att[i-1]
        for j in range(3):
            if att_diff[j] > np.pi:
                att_diff[j] -= 2*np.pi
            if att_diff[j] < -np.pi:
                att_diff[j] += 2*np.pi
        omega[i] = att_diff / dt

    omega[0] = omega[1] - (omega[2] - omega[1])

    for i in range(1, n):
        dt = times[i] - times[i-1]
        domega[i] = (omega[i] - omega[i-1]) / dt

    domega[0] = domega[1] - (domega[2] - domega[1])

    # fit into state vector
    state = np.zeros((n, 20))
    state[:, pos_idx] = positions
    state[:, vel_idx] = v
    state[:, clk_idx] = clock_state
    state[:, att_idx] = att
    state[:, omega_idx] = omega

    # generate imu measurements
    # https://github.com/ethz-asl/kalibr/wiki/IMU-Noise-Model
    acc_bias = np.zeros(3)   # initial bias
    gyro_bias = np.zeros(3)  # initial bias

    imu_meas = np.zeros((n, 6))

    for i in range(1, n):
        # bias as brownian motion
        dt = times[i] - times[i-1]
        acc_bias = acc_bias + imu_model.sigma_ba * np.random.normal(size=3) * np.sqrt(dt)   # random walk sigma_ba = m/s^3 1/sqrt(Hz)
        gyro_bias = gyro_bias + imu_model.sigma_bg * np.random.normal(size=3) * np.sqrt(dt)  # random walk sigma_bg = rad/s^2 1/sqrt(Hz)
        acc_noise = np.random.normal(size=3) * imu_model.sigma_a / np.sqrt(dt)              # sigma_a = m/s^2 1/sqrt(Hz)
        gyro_noise = np.random.normal(size=3) * imu_model.sigma_g / np.sqrt(dt)             # sigma_g = rad/s 1/sqrt(Hz)

        roll, pitch, yaw = att[i]
        R = euler_to_rot(roll, pitch, yaw)
        acc_body = R.T @ a[i]
        acc_meas = acc_body + acc_bias + acc_noise
        gyro_meas = domega[i] + gyro_bias + gyro_noise

        # store bias as states
        state[i, bias_acc_imu_idx] = acc_bias
        state[i, bias_gyro_imu_idx] = gyro_bias
        imu_meas[i] = np.hstack((acc_meas, gyro_meas))

    imu_meas[0] = imu_meas[1]

    return state, imu_meas


def get_gps_measurement(t, state, pos_sat, sigma):

    nsat = pos_sat.shape[0]
    C = 299792458.0  # speed of light

    meas = np.zeros(nsat)

    for i in range(nsat):
        # compute the range
        delta = pos_sat[i] - state[0:3]
        pr = np.linalg.norm(delta) + state[6] * C + np.random.normal() * sigma
        meas[i] = pr

    H = np.zeros((nsat, 20))
    for i in range(nsat):
        delta = pos_sat[i] - state[0:3]
        r = np.linalg.norm(delta)
        H[i, 0:3] = -delta / r
        H[i, 6] = C

    return meas, H