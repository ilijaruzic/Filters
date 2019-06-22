# C filter, full-wave rectifier, comparasion
# All voltages normalized by Vm
# All currents normalizes by wCVm
# Vin = Vm * cos(wt)

from pylab import *

# Initialization
npoint = 720  # reduce if the simulation is too slow
i = arange(npoint)
deg = 360 * (i - 0.5) / npoint
wt = radians(deg)
phi = wt[0 : int(npoint/2)]
minp = cos(wt)
mx = abs(minp)
s = sign(minp)
wjc = -sin(wt) * s

nf = 100  # number of frames (reduce if simulation is too slow)
Joutmax = 0.2

Jout = empty(nf)
Mout = empty(nf)
Moutideal = empty(nf)

close('all')
rc('text', usetex=True)
rc('font', family='serif')
rc('font', size=12)
rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

figure(0, figsize=(4, 4))
degtick = (0, 60, 120, 180, 240, 300, 360)
for k in range(nf):
    jout = (k - 0.5) / nf * Joutmax
    Jout[k] = jout

    beta = arcsin(jout)
    w1 = sqrt(1 - jout**2) - jout * (phi - beta)
    moutideal1 = 1 - jout * phi
    state1 = (phi > beta) & (w1 > mx[0 : int(npoint/2)])
    w = concatenate((w1, w1))
    moutideal = concatenate((moutideal1, moutideal1))
    state = concatenate((state1, state1))

    mout = w * state + mx * ~state
    Mout[k] = mean(mout)
    Moutideal[k] = mean(moutideal)

    clf()
    plot(deg, mout, 'b')
    plot(deg, moutideal, 'r')
    xlabel(r'$\omega t [^\circ]$')
    ylabel(r'$m_{OUT} \, [\text{b}], \; m_{OUT \, ideal} \, [\text{r}]$')
    xlim(0, 360)
    xticks(degtick)
    ylim(0, 1.1)
    savefig('clfilter_comparasion{:d}.png'.format(k), bbox_inches='tight')

figure(1)
plot(Jout, Mout, 'b')
plot(Jout, Moutideal, 'r')
xlabel(r'$J_{OUT}$')
ylabel(r'$M_{OUT} \, [\text{b}], \; M_{OUT \, ideal} \, [\text{r}]$')
xlim(0, Joutmax)
ylim(0, 1.1)
savefig('Mout_comparasion.png', bbox_inches='tight')
