# Potjans & Diesmann cortical microcircuit model implemented on Brian2.

from brian2 import *
import numpy as np
#set_device('genn')

tsim = 1*second

#Number of neurons per layer
#          2/3e   2/3i  4e     4i    5e    5i    6e     6i
n_layer = [20683, 5834, 21915, 5479, 4850, 1065, 14395, 2948]
# Background number per layer
bg_layer = [1600, 1500 ,2100, 1900, 2000, 1900, 2900, 1850]
# Prob. connection table

table = np.array([
		[ 0.101,	0.169,	0.044,	0.082,	0.032,	0.0,	0.008,	0.0 ],
		[ 0.135,	0.137,	0.032,	0.052,	0.075,	0.0,    0.004,	0.0 ],
		[ 0.008,	0.006,	0.050,	0.135,	0.007,	0.0003,	0.045,	0.0 ],
		[ 0.069,	0.003,	0.079,	0.160,	0.003,	0.0,	0.106,	0.0 ],
		[ 0.100,	0.062,	0.051,	0.006,	0.083,	0.373,	0.020,	0.0 ],
		[ 0.055,	0.027,	0.026,	0.002,	0.060,	0.316,	0.009,	0.0 ],
		[ 0.016,	0.007,	0.021,	0.017,	0.057,	0.020,	0.040,	0.225 ],
		[ 0.036,	0.001,	0.003,	0.001,	0.028,	0.008,	0.066,	0.144 ], ])

d_ex = 1.5*ms   # Excitatory delay
d_in = 0.8*ms   # Inhibitory delay

w_ex = 0.15*mV#87.8*pA  # Excitatory weight
g = -4.0        # Inhibitory weight balance

# Neuron model parameters
tau_m   = 10.0*ms   # Membrane time constant
tau_ref = 2.0*ms    # Absolute refractory period
tau_syn = 0.5*ms    # Post-synaptic current time constant
R       = 40*Mohm    # Membrane capacity 
v_r     = -65*mV    # Reset potential
v_th    = -50*mV    # Threshold potential

eqs = '''
	dv/dt = ( -v + v_r + R*I ) / tau_m : volt (unless refractory)
	I : amp
'''

p = [] # Stores NeuronGroups, one for each population

for r in range(0, 8):
	# Addin
	p.append( NeuronGroup(n_layer[r], eqs, threshold='v > v_th', reset = 'v = v_r', method='euler', refractory=tau_ref) )
	p[r].v = v_r # Initialize membrane potencial
	p[r].I = 0.0*pA

# Creatting connections
con = [] # Stores connectiosn
for r in range(0, 8):
	for c in range(0, 8):
		if r%2 == 0:
			w = w_ex
			d = d_ex
		else:
			w = g*w_ex
			d = d_in

		con.append( Synapses(p[r], p[c], on_pre='v_post += w') )
		if r == c:
			con[c+8*r].connect(condition='i!=j', p=table[r][c])
		else:
			con[c+8*r].connect(p=table[r][c])
		con[c+8*r].delay = d


# Creanting BG inputs
bg_in  = []
for r in range(0, 8):
	bg_in.append( PoissonInput(p[r], 'v', bg_layer[r], 8*Hz, weight=w_ex ) )

# Creating Spike Monitors
spikemon = []
for r in range(0, 8):
	spikemon.append( SpikeMonitor( p[r]) )

# Layers frequences
f_23e = spikemon[0].count;
f_23i = spikemon[1].count;
f_4e  = spikemon[2].count;
f_4i  = spikemon[3].count;
f_5e  = spikemon[4].count;
f_5i  = spikemon[5].count;
f_6e  = spikemon[6].count;
f_6i  = spikemon[7].count;

#run(tsim,report='stdout')

print mean(f_23e), mean(f_23i), mean(f_4e), mean(f_4i), mean(f_5e), mean(f_5i), mean(f_6e), mean(f_6i) 

plot(spikemon[0].t/ms, spikemon[0].i,'.k')
show()