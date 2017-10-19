/*----------------------------------------------------------------------------

	Simplified model of bursting cortical neuron
	============================================

        Single-compartment model of "rebound bursts" in pyramidal
        neurons (type of cell very common in association areas of
	cortex).  The model is based on the presence of four
        voltage-dependent currents: 
        - INa, IK: action potentials
        - IM: slow K+ current for spike-frequency adaptation
        - IL: L-type calcium currents for burst generation

  Model described in:

   Pospischil, M., Toledo-Rodriguez, M., Monier, C., Piwkowska, Z., 
   Bal, T., Fregnac, Y., Markram, H. and Destexhe, A.
   Minimal Hodgkin-Huxley type models for different classes of
   cortical and thalamic neurons.
   Biological Cybernetics 99: 427-441, 2008.


        Alain Destexhe, CNRS, 2009
	http://cns.iaf.cnrs-gif.fr

----------------------------------------------------------------------------*/


//----------------------------------------------------------------------------
//  load and define general graphical procedures
//----------------------------------------------------------------------------

load_file("stdrun.hoc")
ismenu = 0
batch = 0
objectvar g[20]			// max 20 graphs
ngraph = 0

proc addgraph() { local ii	// define subroutine to add a new graph
				// addgraph("variable", minvalue, maxvalue)
	ngraph = ngraph+1
	ii = ngraph-1
	g[ii] = new Graph()
	g[ii].size(tstart,tstop,$2,$3)
	//g[ii].size(20000,20050,$2,$3)
	g[ii].xaxis()
	g[ii].yaxis()
	g[ii].addvar($s1,1,0)
	g[ii].save_name("graphList[0].")
	graphList[0].append(g[ii])
}

proc addtext() { local ii	// define subroutine to add a text graph
				// addtext("text")
	ngraph = ngraph+1
	ii = ngraph-1
	g[ii] = new Graph()
	g[ii].size(0,tstop,0,1)
	g[ii].xaxis(3)
	g[ii].yaxis(3)
	g[ii].label(0.1,0.8,$s1)
	g[ii].save_name("graphList[0].")
	graphList[0].append(g[ii])
	text_id = ii
}

proc addline() {		// to add a comment to the text window
				// addline("text")
	g[text_id].label($s1)
}


if(ismenu==0) {
  nrnmainmenu()			// create main menu
  nrncontrolmenu()		// crate control menu
}


//----------------------------------------------------------------------------
//  transient time
//----------------------------------------------------------------------------

trans = 0000

print " "
print ">> Transient time of ",trans," ms"
print " "









//----------------------------------------------------------------------------
//  create PY cells
//----------------------------------------------------------------------------

print " "
print "<<==================================>>"
print "<<            CREATE CELLS          >>"
print "<<==================================>>"
print " "

xopen("TemplatePyramidalCell.hoc")		// read geometry file
xopen("TemplateInhibitoryCell.hoc")

ncells = 3			// nb of cells in each layer <<>>
icells = 1

objectvar IN[icells], PY[ncells], RecurrentSynapse[ncells * ncells-1 + icells], RecurrentConnection[ncells * ncells-1 + icells], InhibitorySynapse[ncells], InhibitoryConnection[ncells] 
for i=0,ncells-1 {
  PY[i] = new sPYb()
}
for i=0,icells-1 {
  IN[i] = new InhibitoryCell()
}

// Synaptic connections //

//pyramidal to pyramidal
x = r.uniform(0,0.02)
synapseCounter = 0
for i = 0, ncells-1{
	access PY[i].soma
	for j = 0, ncells-1{
		if (j!=i){
			PY[i].soma RecurrentSynapse[synapseCounter] = new Exp2Syn()
			RecurrentSynapse[synapseCounter].loc(0.5)
			RecurrentSynapse[synapseCounter].g = 0.0005			//Absolutely no change in output regardless of change in value 
			RecurrentSynapse[synapseCounter].tau1 = 0.5      //rise time of the pulse going through the synapse. originally set to 0.5. see Symes and Wennekers
			RecurrentSynapse[synapseCounter].tau2 = 3 //3		//fall time of the pulse going through the synapse. originally set to 3
			RecurrentSynapse[synapseCounter].e = 0
			PY[j].soma RecurrentConnection[synapseCounter] = new NetCon(&PY[j].soma.v(.5), RecurrentSynapse[synapseCounter], 0, 2, x)
			synapseCounter = synapseCounter + 1
		}	
	}
}

// inhibitory to pyramidal
gabacounter = 0
for i = 0, ncells-1{
	access PY[i].soma
	for a = 0, icells -1{
		
		PY[i].soma InhibitorySynapse[gabacounter] = new Exp2Syn() 
		InhibitorySynapse[gabacounter].loc(0.5)
		InhibitorySynapse[gabacounter].g = 1			
		InhibitorySynapse[gabacounter].tau1 = 1      
		InhibitorySynapse[gabacounter].tau2 = 250		
		InhibitorySynapse[gabacounter].e = -75
		//IN[a].soma InhibitoryConnection[gabacounter] = new NetCon(&IN[a].soma.v(.5), InhibitorySynapse[gabacounter], 0, 2, 0.02)
		IN[a].soma InhibitoryConnection[gabacounter] = new NetCon(&IN[a].soma.v(.5), InhibitorySynapse[gabacounter], 0, 2, 0) //0.01111)
		
		gabacounter = gabacounter + 1	
	}
}

// pyramidal to inhibitory
for i = 0, ncells-1{
	for a = 0, icells -1{
		access IN[a].soma
		IN[a].soma RecurrentSynapse[synapseCounter] = new Exp2Syn() //neuron to interneuron
		RecurrentSynapse[synapseCounter].loc(0.5)
		RecurrentSynapse[synapseCounter].g = 0.0005			//Absolutely no change in output regardless of change in value 
		RecurrentSynapse[synapseCounter].tau1 = 0.5      //rise time of the pulse going through the synapse. originally set to 0.5. see Symes and Wennekers
		RecurrentSynapse[synapseCounter].tau2 = 3		//fall time of the pulse going through the synapse. originally set to 3
		RecurrentSynapse[synapseCounter].e = 0
		PY[i].soma RecurrentConnection[synapseCounter] = new NetCon(&PY[i].soma.v(.5), RecurrentSynapse[synapseCounter], 0, 2, 0) //0.02)
		synapseCounter = synapseCounter + 1
		print gabacounter

	}
}






//----------------------------------------------------------------------------
//  insert electrode in each PY cell
//----------------------------------------------------------------------------

if(ismenu==0) {
  load_file("electrod.hoc")	// electrode template
  ismenu = 1
}

objectvar El[2 * ncells]			// create electrodes

CURR_AMP = 0.15

//for i=0,ncells-1 {			// Positive pulse current
for i=0,1 {			// Positive pulse current
	PY[i].soma El[i] = new Electrode()
	PY[i].soma El[i].stim.loc(0.5)
	El[i].stim.del = 5000
	El[i].stim.dur = 2000
	El[i].stim.amp = CURR_AMP
}

for i=0,ncells-1 {			// Negative pulse current
	PY[i].soma El[i +ncells] = new Electrode()
	PY[i].soma El[i +ncells].stim.loc(0.5)
	El[i+ncells].stim.del = 12000
	El[i+ncells].stim.dur = 0  //15
	El[i+ncells].stim.amp = -0.0
}

objref netstim[ncells], RandSynapse[ncells], RandConnections[ncells], r
for i = 0, ncells-1{         // Sine wave injection
	access PY[i].soma
	PY[i].soma netstim[i] = new SinClamp(.5)
	netstim[i].del = 12000
	netstim[i].dur = 280000
	netstim[i].freq = 7
	netstim[i].pkamp = 0.0
}

//for i=0,icells-1 {			// Stimulation to inhibitory cells
	//IN[i].soma El[ncells + i] = new Electrode()
	//IN[i].soma El[ncells + i].stim.loc(0.5)
	//El[ncells + i].stim.del = 3000
	//El[ncells + i].stim.dur = 10
	//El[ncells + i].stim.amp = -0.5
//}

electrodes_present=1



//----------------------------------------------------------------------------
//  setup simulation parameters
//----------------------------------------------------------------------------

Dt = .1				// macroscopic time step <<>>
npoints = 380001    //280000 in Fig4,   380001 in Fig5

dt = 0.1			// must be submultiple of Dt
tstart = trans
tstop = trans + npoints * Dt
runStopAt = tstop
steps_per_ms = 5
celsius = 36
v_init = -84






//----------------------------------------------------------------------------
//  add graphs
//----------------------------------------------------------------------------

strdef gtxt

if(batch == 0) {
  //addgraph("PY[0].soma.cai(0.5)",0,0.016)
  //addgraph("PY[0].soma.gcan_ican(0.5)",0,1.2e-7)
  
//  addgraph("PY[0].soma.m_iahp",0,1)
  for i=0,ncells-1 {
	sprint(gtxt,"PY[%d].soma.v(0.5)",i)
	addgraph(gtxt,-120,40)
  }
}





//----------------------------------------------------------------------------
//  add text
//----------------------------------------------------------------------------

objref ActionPotentialCounter, SpikeVector
SpikeVector = new Vector()
access PY[0].soma
PY[0].soma ActionPotentialCounter = new APCount()
ActionPotentialCounter.loc(0.5)
ActionPotentialCounter.thresh = 0
ActionPotentialCounter.record(SpikeVector)

objref ActionPotentialCounter2, SpikeVector2
SpikeVector2 = new Vector()
access PY[1].soma
PY[1].soma ActionPotentialCounter2 = new APCount()
ActionPotentialCounter2.loc(0.5)
ActionPotentialCounter2.thresh = 0
ActionPotentialCounter2.record(SpikeVector2)

objref ActionPotentialCounter3, SpikeVector3
SpikeVector3 = new Vector()
access PY[2].soma
PY[2].soma ActionPotentialCounter3 = new APCount()
ActionPotentialCounter3.loc(0.5)
ActionPotentialCounter3.thresh = 0
ActionPotentialCounter3.record(SpikeVector3)

min_p1 = 1.0   // Parameter 1 range setting (gcan)
delta_p1 = 1  // increment
stepCount_p1 = 1  // number of steps

min_p2 = 0.0   // Parameter 2 range setting (Wpp)
delta_p2 = 0.001
stepCount_p2 = 1

objref DeltaVector, DeltaVector2, ResultVector, ResultVector2, resultGraph, outputFile1, outputFile2, outputFile3

DeltaVector = new Vector(stepCount_p1)
DeltaVector2 = new Vector(stepCount_p2)


proc RunExperiment() {
	
	outputFile1 = new File()
	outputFile1.wopen("Result_freq.dat")
	
	// Can current on off
	for i=0,ncells-1 {
		PY[i].soma.gbar_ican = 8.67e-6 * 1
	}
	
	for stepCounter = 0, stepCount_p1-1 {
		p1 = min_p1 + (stepCounter * delta_p1)
		DeltaVector.x[stepCounter] = p1
		
		// Parameter 1 /////////////////////////////////////
		//for c = 0, (ncells * (ncells-1))-1 {
			//RecurrentSynapse[c].tau2 = p1
		//}
		
		for i=0,ncells-1 {
			PY[i].soma.gbar_ican = 8.67e-6 * p1
		}
		
		//for i=0,ncells-1 {
				//El[i+ncells].stim.amp = p1
		//}
		
		//for i=0,ncells-1 {			
		//		InhibitorySynapse[i].tau2 = p1
		//}
		
		for stepCounter2 = 0 , stepCount_p2-1 {
			p2 = min_p2 + (stepCounter2 * delta_p2)
			DeltaVector2.x[stepCounter2] = p2
			
			// Parameter 2 /////////////////////////////////////
			for i=0,(ncells * (ncells-1))-1 {	// This is corrrect
				RecurrentConnection[i].weight = p2
			}
			
			//for i=0,ncells-1 {			
				//El[i+ncells].stim.dur = p2
			//}
			
			//for i = 0, ncells-1 {
			//	netstim[i].pkamp = p2
			//}
			
			//for i=0,ncells-1 {
				//InhibitoryConnection[i].weight[0] = p2
			//}
		
			run()
			
			currentFreq = 0	
			currentFreq2 = 0
			currentFreq3 = 0
			
			// Normal cases
			//totalDelay = El[0].stim.del + El[0].stim.dur + 10000
			//totalInterval = El[0].stim.del + El[0].stim.dur + 20000
			
			// In case of netative current injection
			totalDelay = El[ncells + 1].stim.del + El[ncells + 1].stim.dur + 5000
			totalInterval = totalDelay + 10000
			
				for q = 0, SpikeVector.size-1 {
					currentValue = SpikeVector.x[q]
					if (currentValue > totalDelay){
						if (SpikeVector.x[q] < totalInterval){
							currentFreq = currentFreq + 1
						}
					}
				}
				for q = 0, SpikeVector2.size-1 {
					currentValue = SpikeVector2.x[q]
					if (currentValue > totalDelay){
						if (SpikeVector2.x[q] < totalInterval){
							currentFreq2 = currentFreq2 + 1
						}
					}
				}
				for q = 0, SpikeVector3.size-1 {
				currentValue = SpikeVector3.x[q]
				if (currentValue > totalDelay){
					if (SpikeVector3.x[q] < totalInterval){
						currentFreq3 = currentFreq3 + 1
					}
				}
			}
						
			currentFreq = currentFreq/10
			currentFreq2 = currentFreq2/10
			currentFreq3 = currentFreq3/10
			
			progress = ((stepCounter * stepCount_p2) + (stepCounter2 + 1))/(stepCount_p1 * stepCount_p2) * 100
			print " "
			print "******************************************"
			print "Progress in percentage"
			print progress
			print "p1"
			print p1
			print "p2"
			print p2
			print "First cell freq"
			print currentFreq
			print "Second cell freq"
			print currentFreq2
			print "Third cell freq"
			print currentFreq3
			
			outputFile1.printf("%f", currentFreq2)
			outputFile1.printf(",")

		}
		outputFile1.printf("\n")

	}
	outputFile2 = new File()
	outputFile2.wopen("Result_x.dat")
	for j = 0, stepCount_p1-1 {
		outputFile2.printf("%f", DeltaVector.x[j])
		outputFile2.printf("\n")
	}
	outputFile3 = new File()
	outputFile3.wopen("Result_y.dat")
	for j = 0, stepCount_p2-1 {
		outputFile3.printf("%f", DeltaVector2.x[j])
		outputFile3.printf("\n")
	}
	
	outputFile1.close()
	outputFile2.close()
	outputFile3.close()
	print "--------   Finished simulation"
	
}


isICAN = 1


proc ShowButtonSet() {
  xpanel("Test panel")
  xvalue("Calcium decay", "PY[0].soma.taur_cad(0.5)")
  xbutton("Run Complete experiment", "RunExperiment()") 
  xvalue("gbar Im", "PY[0].soma.gkbar_im(0.5)")
  xvalue("Recurrent synapse weight", "RecurrentConnection.weight")  
  xcheckbox("Ican", &isICAN,  "SetIcan()")
  
  xpanel(100, 800)
}

proc SetIcan(){
	if (isICAN == 0){		
		for i=0,ncells-1 {
			access PY[i].soma
			uninsert ican
		}
	}
	if (isICAN == 1){
		for i=0,ncells-1 {
			access PY[i].soma
			insert ican
		}
	}
}

load_file("data_saving.hoc")
ShowButtonSet()
SetIcan()
