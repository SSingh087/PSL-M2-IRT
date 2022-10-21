import numpy as np
import matplotlib.pyplot as plt

xmin, xmax = 0, 1
ymin, ymax = 0, 1
tmin, tmax = 0, 5

# Polytropic index 
gamma = 5/3

x_all_cells, y_all_cells, t_all_cells = 200, 10, 20
dx = (xmax - xmin) / x_all_cells
delta_y = (ymax - ymin) / x_all_cells

def time_step_get(gamma, rho, vx, vy, p, dx):
    cmax = characteristic_speed_get(gamma, rho, vx, vy, p)
    dt = dx/np.amax((cmax))
    return dt

def characteristic_speed_get(gamma, rho, v, p):
    # v can be vx or vy
    csound = csound_get(gamma, rho, p)
    cmax = v + csound
    cmin = v - csound
    return cmin, cmax

def characteristic_speed_norm_get(gamma,rho, vx, vy, p):
    csound = csound_get(gamma, rho, p)
    cmax = csound + np.sqrt(vx**2 + vy**2)
    return cmax

def csound_get(gamma, rho, p):
    return np.sqrt(gamma * p / rho)

def minmod(a, b):
    1/2*(np.sign(a)+np.sign(b))*min(abs(a), abs(b))

def minmod_build(var, varL ,varR):
    var_intL = var + 0.5 * minmod(varR - var, var - varL)
    var_intR = varR - 0.5 * minmod(np.roll(varR, 1) - varR, varR -var)
    return var_intL,var_intR

def space_reconstrution(idir, rho, vx, vy, p):
    i_shift = 1
    rho_R = np.roll(rho,i_shift)
    vx_R = np.roll(vx,i_shift)
    vy_R = np.roll(vy,i_shift)
    p_R = np.roll(p,i_shift)
    i_shift = -1
    rho_L = np.roll(rho, i_shift)
    vx_L = np.roll(vx, i_shift)
    vy_L = np.roll(vy, i_shift)
    p_L = np.roll(p, i_shift)
    rho_intL, rho_intR = minmod_build(rho, rho_L, rho_R)
    vx_intL, vx_intR = minmod_build(vx, vx_L, vx_R)
    vy_intL, vy_intR = minmod_build(vy, vy_L, vy_R)
    p_intL, p_intR = minmod_build(p, p_L, p_R)
    return rho_intL, rho_intR, vx_intL, vx_intR, vy_intL, vy_intR, \
            p_intL, p_intR


def flux_get(idir, gamma, rho, vx, vy, p):
    if idir == 1 :
        flux_rho = rho*vx
        flux_momx=p + rho*vx**2
        flux_momy=rho * vx **2
        flux_e = (p / (gamma - 1) + rho * np.sqrt(vx**2 + vy**2) / 2.0) * vx
    else :
        flux_rho = rho * vy
        flux_momx = rho * vx**2
        flux_momy = p+rho*vy**2
        flux_e = (p / (gamma - 1) + rho * np.sqrt(vx**2 + vy**2) / 2.0) * vy
        return flux_rho,flux_momx,flux_momy,flux_e

def solver_hll_get_flux(varL, varR, sL, sR, fluxL, fluxR):
    if sL >0:
        flux = fluxL
    elif sR < 0:
        flux=fluxR
    else :
        flux=(sR * fluxL - sL * fluxR + sL * sR * (varR - varL)) / (sR - sL)
    return flux

def advance(dt, gamma, dx, U_int, U_old):
    Dflux_total = np.zeros((linX,linY))
    for idir in [0,1]: # loop direction
        U_int_prim = phys_prim(gamma, U_int[0], U_int[3], U_int[4], U_int[5])
        U_int_prim_surfL, U_int_prim_surfR = space_reconstrution(i, 
             U_int_prim[0], U_int_prim[1], U_int_prim[2], U_int_prim[6])
        sL = characteristic_speed_get(gamma,
             U_int_prim_surfL[0], U_int_prim_surfL[1], U_int_prim_surfL[2], U_int_prim_surfL[6])
        sR = characteristic_speed_get(gamma, 
             U_int_prim_surfR[0], U_int_prim_surfR[1], U_int_prim_surfR[2], U_int_prim_surfR[6])
        flux_surfL = flux_get(idir, gamma,
             U_int_prim_surfL[0], U_int_prim_surfL[1], U_int_prim_surfL[2], U_int_prim_surfL[6])
        flux_surfR = flux_get(idir, gamma,
             U_int_prim_surfR[0], U_int_prim_surfR[1], U_int_prim_surfR[2], U_int_prim_surfR[6])
        flux_hll = solver_hll_get_flux(U_int_prim_surfL,
               U_int_prim_surfR, sL, sR, flux_surfL, flux_surfR)
        Dflux_total += (flux_hll - np.roll(flux_hll, -1, axis=idir))/dx
    U = U_old - dt * Dflux_total
    return U

def phys_prim(gamma, rho, momx, momy, energy):
    vx = momx / rho
    vy = momy / rho
    p = (energy - 0.5 * rho * (vx**2 + vy**2)) * (gamma - 1)
    return rho, vx, vy, p

def phys_cons(gamma, rho, vx, vy, p):
    momx = rho * vx
    momy = rho * vy
    energy = p / (gamma - 1.) + 0.5 * rho
    return rho, momx, momy, energy

if __name__ == "__main__":
    
    T = np.linspace(tmin, tmax, t_all_cells, endpoint=True)
    linX = np.linspace(xmin-dx/2, xmax+dx/2, x_all_cells, endpoint=True)
    linY = np.linspace(ymin-dx/2, ymax+dx/2, y_all_cells, endpoint=True)
    x, y = np.meshgrid(linX, linY)
    
    # rho
    rho = np.ones(x.shape)
    
    # pressure
    p = 0.1*np.ones(x.shape)
    
    # velocity
    vx, vy = np.ones(x.shape), np.ones(x.shape)

    # declare initial velocity
    for j in range(len(linY)):
        for i in range(len(linX)):
            if linY[j] < (ymax-ymin)/2:
                vx[j, i] = -0.5
            elif linY[j] > (ymax-ymin)/2:
                vx[j, i] = 0.5 
            vy[j, i] = -0.01*np.sin(2*np.pi*((linX[i]-xmin)/(xmax-xmin)))* \
                (np.exp(-(np.abs(linY[j]-(ymax-ymin)/2))/(2. * 0.1 )))
    
    for t in range(len(T)):
        # evaluate conservative variables.
        rho, momx, momy, energy = phys_cons(gamma, rho, vx, vy, p)

        # copy conservative variables to old variables
        rho_old = np.copy(rho)
        momx_old = np.copy(momx)
        momy_old = np.copy(momy)
        e_old = np.copy(energy)
        
        # evaluate primitive variables from conservative variables
        rho, vx, vy, p = phys_prim(gamma, rho_old, momx_old, momy_old, e_old)
        
        # Compute the next timestep Δt
        dt = time_step_get(gamma, rho, vx, vy, p, dx)
        U = np.array([rho, vx, vy, momx, momy, energy, p])
        U_old = np.array([rho_old, vx, vy, momx_old, momy_old, e_old, p])
        U_int = np.copy(U)
        U = advance(dt/2, gamma, dx, U_int, U_old)
        U = advance(dt, gamma, dx, U, U_old)
        for idir in [0,1]:
            rhoL, rhoR, vxL, vxR, vyL, vyR, pL, pR = space_reconstrution(idir, )
            cmin_left, cmax_left = characteristic_speed_get(gamma, rhoL, vL, pL)
        cmin_right, cmax_right = characteristic_speed_get(gamma, rhoR, vR, pR)
        cmin=min(cmin_left,cmin_right)
        cmax=max(cmax_right,cmax_right)
        
