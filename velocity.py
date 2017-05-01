from params import *


def termalAccomodation(speed, bound):
    """Calculate new temperature based on
     accomodation coefficient"""
    if bound in ['up', 'right']:
        T_wall = T_w
    else:
        T_wall = T
    T_i = m * speed ** 2 / 3 / k
    T_r = T_i * (1 - alpha) + alpha * T_wall
    beta_new = sqrt(m / 2 / k / T_r)
    return beta_new


def acceptance_rejection(s, beta):
    """Acceptance-rejection procedure"""

    def f_bu(u, s, beta):
        return (beta * u + s) * np.exp(-beta ** 2 * u ** 2)

    f_max = (s + sqrt(s ** 2 + 2)) / 2 / exp(0.5 + s / 2 * (s - sqrt(s ** 2 + 2)))
    a, b = 0, 3 / beta  # cut-off for pdf
    while True:
        R_f = rand.random()
        u = rand.uniform(a, b)
        if R_f < f_bu(u, s, beta) / f_max:
            return u


def diffuseVelocity(beta_new, bound):
    """Generate velocities for diffuse collision
    using acceptance-rejection method for component
    normal to surface and direct pairs for others"""
    phi = rand.uniform(0, 2 * pi)
    rho = sqrt(-log(rand.random())) / beta_new
    if bound in ['left', 'right']:
        u_y = rho * cos(phi)
        u_z = rho * sin(phi)
        if bound == 'left':
            u_x = acceptance_rejection(0, beta_new)
        else:
            u_x = -acceptance_rejection(0, beta_new)
    if bound in ['down', 'up']:
        u_x = rho * cos(phi)
        u_z = rho * sin(phi)
        if bound == 'down':
            u_y = acceptance_rejection(0, beta_new)
        else:
            u_y = -acceptance_rejection(0, beta_new)
    return [u_x, u_y, u_z]


def newVelocity(bound, U):
    """Calculate reflected velocities
    after collision with the box wall"""
    if reflection_type == 'specular':  # switch normal component
        if bound in ['down', 'up']:
            U[1] = -U[1]
        elif bound in ['left', 'right']:
            U[0] = -U[0]
        else:
            sys.exit("Molecule:newVelocity(): wrong bound")
    elif reflection_type == 'diffuse':
        speed = norm(U)
        beta_new = termalAccomodation(speed, bound)
        U = diffuseVelocity(beta_new, bound)
    if COUETTE and (bound is 'up'):
        U[0] += u_couette
    return U[0], U[1], U[2]
