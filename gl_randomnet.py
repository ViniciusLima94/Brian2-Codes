# A simple random network with GL model

from brian2 import *
import matplotlib 

defaultclock.dt = 1.0*ms

n_ex = 800  # Number of excitatory neurons
n_in = 200  # Number of inhibitory neurons
n = n_ex + n_in # Total number of neurons

tsim = 5.0*second # simulation time

we = 0.5*mV;  # Excitatory weigth
g = -4.0;  # Balance for inhibitory weigth

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

# GL_RS Parameters
R_rs = 10e-3*Gohm  # Resistence
tau_rs = 10*ms     # Time constant
vt_rs = -64.5*mV   # Threshold potential
vs_rs = -37.5*mV   # Saturation potential
# phi_v parameters
gamma_rs = 0.037/mV 
r_rs = 1.0
delta_rs = 0.0

# GL_FS Parameters
R_fs = 16e-3*Gohm  # Resistence
tau_fs = 10*ms     # Time constant
vt_fs = -63.1*mV   # Threshold potential
vs_fs = -50.50*mV   # Saturation potential
# phi_v parameters
gamma_fs = 0.05/mV 
r_fs = 1.0
delta_fs = 0.3

vr = -65*mV     # Reset potential

# Excitatory population
p_ex = NeuronGroup(n_ex, 
''' 
dv_rs/dt = ( -v_rs + vr + R_rs*I ) / tau_rs : volt (unless refractory)
I : amp
''', 
	              threshold='rand() < phi_v(v_rs,vt_rs, vs_rs, vr, gamma_rs,r_rs,delta_rs)', 
	              reset='v_rs=vr', 
	              method='euler', 
	              refractory=1.0*ms)

p_ex.v_rs = 'rand()*(vs_rs-vt_rs) + vt_rs'
p_ex.I = 30*pA

# Inhibitory population
p_in = NeuronGroup(n_in, 
''' 
dv_fs/dt = ( -v_fs + vr + R_fs*I ) / tau_fs : volt (unless refractory)
I : amp
''', 
	              threshold='rand() < phi_v(v_fs,vt_fs, vs_fs, vr, gamma_fs,r_fs,delta_fs)', 
	              reset='v_fs=vr', 
	              method='euler', 
	              refractory=1.0*ms)

p_in.v_fs = 'rand()*(vs_fs-vt_fs) + vt_fs'
p_in.I = 30*pA

# Connections EX-> EX
con_ex_ex = Synapses(p_ex, p_ex, on_pre='v_rs_post += we')
con_ex_ex.connect(condition='i!=j', p=0.1)
con_ex_ex.delay = 1.5*ms

# Connections EX -> IN
con_ex_in = Synapses(p_ex, p_in, on_pre='v_fs_post += we')
con_ex_in.connect(p=0.1)
con_ex_in.delay = 1.5*ms

# Connections IN -> EX
con_in_ex = Synapses(p_in, p_ex, on_pre='v_rs_post += g*we')
con_in_ex.connect(p=0.1)
con_in_ex.delay = 1.5*ms

# Connections IN -> IN
con_in_in = Synapses(p_in, p_in, on_pre='v_fs_post += g*we')
con_in_in.connect(condition='i!=j', p=0.1)
con_in_in.delay = 1.5*ms

spikemon1 = SpikeMonitor(p_ex)
spikemon2 = SpikeMonitor(p_in)

freq1 = spikemon1.count;
freq2 = spikemon2.count;

run(tsim, report='stdout')

print mean(freq1), mean(freq2)


plot(spikemon1.t/ms, spikemon1.i,'.k')
plot(spikemon2.t/ms, spikemon2.i+n_ex,'.r')
xlabel('Time (ms)')
ylabel('Neuron index');
show()