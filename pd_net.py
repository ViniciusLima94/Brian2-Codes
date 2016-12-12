# Potjans & Diesmann cortical microcircuit model implemented on Brian2.

from brian2 import *
import numpy as np
#set_device('genn')

import matplotlib
matplotlib.use('Agg')

set_device('cpp_standalone', directory='PD')
prefs.devices.cpp_standalone.openmp_threads = 8

defaultclock.dt = 0.1*ms

tsim = 5*second

N = 80000 # Total population
# Fraction of neurons in each layer
#       2/3e   2/3i   4e     4i     5e     5i     6e     6i
frac = [.2680, .0756, .2840, .0710, .0628, .0138, .1866, .0382 ]

#Number of neurons per layer
n_layer = [N*f for f in frac]
n_layer = [int(round(n_pop)) for n_pop in n_layer]
# Reescale factor
rf = 80000.0 / N 
nn_cum = cumsum(n_layer)
# Background number per layer
bg_layer = [1600, 1500 ,2100, 1900, 2000, 1900, 2900, 1850]
bg_layer = [bg_pop/rf for bg_pop in bg_layer]   # Reescaling
# Prob. connection table

table = [[0.1009,  0.1689, 0.0437, 0.0818, 0.0323, 0.,     0.0076, 0.    ],
        [0.1346,   0.1371, 0.0316, 0.0515, 0.0755, 0.,     0.0042, 0.    ],
        [0.0077,   0.0059, 0.0497, 0.135,  0.0067, 0.0003, 0.0453, 0.    ],
        [0.0691,   0.0029, 0.0794, 0.1597, 0.0033, 0.,     0.1057, 0.    ],
        [0.1004,   0.0622, 0.0505, 0.0057, 0.0831, 0.3726, 0.0204, 0.    ],
        [0.0548,   0.0269, 0.0257, 0.0022, 0.06,   0.3158, 0.0086, 0.    ],
        [0.0156,   0.0066, 0.0211, 0.0166, 0.0572, 0.0197, 0.0396, 0.2252],
        [0.0364,   0.001,  0.0034, 0.0005, 0.0277, 0.008,  0.0658, 0.144 ]]

d_ex = 1.5*ms      # Excitatory delay
std_d_ex = 0.75*ms # Std. Excitatory delay
d_in = 0.8*ms      # Inhibitory delay
std_d_in = 0.4*ms  # Std. Inhibitory delay

w_ex = rf*.3512*mV      # Excitatory weight
std_w_ex = rf*.0352*mV  # Standar deviation weigth
g = 4.0                 # Inhibitory weight balance

# Neuron model parameters
tau_m   = 10.0*ms   # Membrane time constant
tau_ref = 2.0*ms    # Absolute refractory period
tau_syn = 0.5*ms    # Post-synaptic current time constant
R       = 40*Mohm   # Membrane capacity 
v_r     = -65*mV    # Reset potential
v_th    = -50*mV    # Threshold potential

eqs = '''
	dv/dt = ( -v + v_r + R*I ) / tau_m : volt (unless refractory)
	I : amp
'''

neurons = NeuronGroup(sum(n_layer), eqs, threshold='v>v_th', reset='v=v_r', method='euler', refractory=tau_ref)
neurons.v = 'rand()*(v_th-v_r) + v_r'
neurons.I = 0.0*pA

p = [] # Stores NeuronGroups, one for each population

for r in range(0, 8):
	if r == 0:
		p.append(neurons[:nn_cum[0]])
	else:
		p.append(neurons[nn_cum[r-1]:nn_cum[r]])

# Creatting connections
con = [] # Stores connections
for c in range(0, 8):
	for r in range(0, 8):
		# Excitatory layer
		if (c % 2) == 0:	
			if c == 2 and r == 0:
				con.append( Synapses(p[c],p[r], """ w:volt """,on_pre = 'v_post += 2*w') )
				con[-1].connect(condition='i!=j', p=table[r][c])
				con[-1].w = 'abs(w_ex + std_w_ex*randn())'
			else:
				con.append( Synapses(p[c],p[r], """ w:volt """,on_pre = 'v_post += w') )
				con[-1].connect(condition='i!=j', p=table[r][c])
				con[-1].w = 'abs(w_ex + std_w_ex*randn())'
			con[-1].delay = 'abs(d_ex + std_d_ex*randn())'
		# Inhibitory layer		
		else:
			con.append( Synapses(p[c],p[r], """ w:volt """,on_pre = 'v_post -= g*w') )
			con[-1].connect(condition='i!=j', p=table[r][c])
			con[-1].w = 'abs(w_ex + std_w_ex*randn())'
           		con[-1].delay = 'abs(d_in + std_d_in*randn())'
			
# Creanting BG inputs
bg_in  = []
for r in range(0, 8):
	bg_in.append( PoissonInput(p[r], 'v', bg_layer[r], 8*Hz, weight = w_ex ) )

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

smon_net = SpikeMonitor(neurons)

net = Network(collect())
net.add(neurons,p, con, spikemon, bg_in)
net.run(tsim,report='stdout')

# Saving some values in files
# Frequences 
f1 = open('fperlayer.dat', 'w')
f1.write(str(mean(f_23e)))
f1.write("\n")
f1.write(str(mean(f_23i)))
f1.write("\n")
f1.write(str(mean(f_4e)))
f1.write("\n")
f1.write(str(mean(f_4i)))
f1.write("\n")
f1.write(str(mean(f_5e)))
f1.write("\n")
f1.write(str(mean(f_5i)))
f1.write("\n")
f1.write(str(mean(f_6e)))
f1.write("\n")
f1.write(str(mean(f_6i)))

# Saving some excitatory weights
f2 = open('wex.dat', 'w')
for wex in con[0].w[:5000]:
	f2.write(str(wex/mV))
	f2.write("\n")

# Saving some inhibitory weights
f3 = open('win.dat', 'w')
for win in con[-1].w[:5000]:
	f3.write(str(win/mV))
	f3.write("\n")

# Saving some excitatory delays
f4 = open('dex.dat', 'w')
for dex in con[0].delay[:5000]:
	f4.write(str(dex/ms))
	f4.write("\n")

# Saving some inhibitory delays
f5 = open('din.dat', 'w')
for din in con[-1].delay[:5000]:
	f5.write(str(din/ms))
	f5.write("\n")

# Saving raster plot
f6 = open('raster.dat', 'w')
for idx in range(0, len(smon_net.t)):
	f6.write(str(smon_net.t[idx]/ms))
	f6.write("\t")
	f6.write(str(smon_net.i[idx]))
	f6.write("\n")

f1.close()
f2.close()
f3.close()
f4.close()
f5.close()
f6.close()

# plot(smon_net.t/ms, smon_net.i,'.k')
# xlabel('Time (ms)')
# ylabel('Neuron index');
# ylim(0,sum(n_layer))
# xlim(500,1000)
# plt.gca().invert_yaxis()
# savefig('raster_PDnet.png')
# show()
