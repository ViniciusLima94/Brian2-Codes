#GL MODEL

from brian2 import *
import matplotlib 

defaultclock.dt = 1.0*ms

@implementation('cpp', '''
     #include <math.h>
     double phi_v(double v, double vt, double vs, double vr, double gamma, double r, double delta) {
        if(v <= vt) {
        	return 0;
        } else if(v > vt && v <= vs) {
        	return pow(gamma*(v-vt), r) + delta;
        } else if(v > vs) {
        	return 1;
        }
     }
     ''')
@check_units(v=volt, vt=volt, vs=volt, vr=volt, gamma=1/volt,r=1,delta=1, result=1)
def phi_v(v,vt, vs, vr, gamma,r,delta):
	if v <= vt:
		return 0
	elif v > vt and v<=vs:
		return (gamma*(v-vt))**r + delta
	elif v > vs:
		return 1

el_class = 'FS';  # Electrophysiologic class (RS or FS)

if el_class == 'RS':
	R = 10e-3*Gohm  # Resistence
	tau = 10*ms     # Time constant
	vr = -65*mV     # Reset potential
	vt = -64.5*mV   # Threshold potential
	vs = -37.5*mV   # Saturation potential
	# phi_v parameters
	gamma = 0.037/mV 
	r = 1.0
	delta = 0.0
elif el_class == 'FS':
	R = 16e-3*Gohm  # Resistence
	tau = 10*ms     # Time constant
	vr = -65*mV     # Reset potential
	vt = -63.1*mV   # Threshold potential
	vs = -50.50*mV   # Saturation potential
	# phi_v parameters
	gamma = 0.05/mV 
	r = 1.0
	delta = 0.3

# Modeled as LIF where mu = (1 - dt/tau), mu = 0.9, dt = 1.0 ms and tau = 10 ms
eqs = '''
	dv/dt = ( -v + vr + R*I ) / tau : volt (unless refractory)
	I : amp
'''

n = NeuronGroup(1, eqs, threshold='rand()< phi_v(v,vt, vs, vr, gamma,r,delta)', reset='v=vr', method='euler', refractory=0.0*ms)
n.v = vr

trace = StateMonitor(n, 'v', record=0)

n.I = 120*pA
run(1000 * ms)

plot(trace.t, trace.v[0] / mV)
xlabel('time (ms)')
ylabel('membrane potential (mV)')
show()
