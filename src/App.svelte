<script>
  
  import { scaleLinear } from "d3-scale";
  // import { Date } from "d3-time"
  import Chart from './Chart.svelte';
  import { onMount } from 'svelte';
  import { selectAll } from 'd3-selection'
  import { drag } from 'd3-drag';
  import queryString from "query-string";
  import Checkbox from './Checkbox.svelte';
  import Arrow from './Arrow.svelte';
  import { format } from 'd3-format'
  import { event } from 'd3-selection'

  import katex from 'katex';

  const legendheight = 67 

  function range(n){
    return Array(n).fill().map((_, i) => i);
  }

  function formatNumber(num) {
    return num.toString().replace(/(\d)(?=(\d{3})+(?!\d))/g, '$1,')
  }

  var sum = function(arr, bools){
    var x = 0
    for (var i = 0; i < arr.length; i++) {
      x = x + arr[i]*(bools[i] ? 1 : 0)
    }
    return x
  }

  var Integrators = {
    Euler    : [[1]],
    Midpoint : [[.5,.5],[0, 1]],
    Heun     : [[1, 1],[.5,.5]],
    Ralston  : [[2/3,2/3],[.25,.75]],
    K3       : [[.5,.5],[1,-1,2],[1/6,2/3,1/6]],
    SSP33    : [[1,1],[.5,.25,.25],[1/6,1/6,2/3]],
    SSP43    : [[.5,.5],[1,.5,.5],[.5,1/6,1/6,1/6],[1/6,1/6,1/6,1/2]],
    RK4      : [[.5,.5],[.5,0,.5],[1,0,0,1],[1/6,1/3,1/3,1/6]],
    RK38     : [[1/3,1/3],[2/3,-1/3,1],[1,1,-1,1],[1/8,3/8,3/8,1/8]]
  };

  // f is a func of time t and state y
  // y is the initial state, t is the time, h is the timestep
  // updated y is returned.
  var integrate=(m,f,y,t,h)=>{
    for (var k=[],ki=0; ki<m.length; ki++) {
      var _y=y.slice(), dt=ki?((m[ki-1][0])*h):0;
      for (var l=0; l<_y.length; l++) for (var j=1; j<=ki; j++) _y[l]=_y[l]+h*(m[ki-1][j])*(k[ki-1][l]);
      k[ki]=f(t+dt,_y,dt); 
    }
    for (var r=y.slice(),l=0; l<_y.length; l++) for (var j=0; j<k.length; j++) r[l]=r[l]+h*(k[j][l])*(m[ki-1][j]);
    return r;
  }

  $: country           = "Belgium"

  $: logN              = Math.log(11515790)
  $: N                 = Math.exp(logN)
  $: I0                = 31
  $: R0                = 3.3
  $: InterventionTime  = 100  
  $: OMInterventionAmt = 0.7
  $: InterventionAmt   = 1 - OMInterventionAmt
  $: Time              = 220
  $: Xmax              = 110000
  $: dt                = 4
  
  $: duration          = 7*12*1e10  // deprecate
  
  $: N_ICU             = 1950 * 1e5 / N 
  
  $: P_symptoms        = 0.6
  $: P_hospital        = 0.15
  $: P_ICU             = 1/6.5
  $: P_die_in_ICU      = 0.4
  $: P_die_in_hospital = 0.0

  $: D_to_asym         = 5
  $: D_to_symptoms     = 2
  $: D_to_hospital     = 7
  $: D_to_ICU          = 3.5
  $: D_to_death        = 4 / Math.log(2) // median

  $: D_recov_asym      = 14
  $: D_recov_symptoms  = 14
  $: D_recov_ICU       = 21 //3/ Math.log(2) // median
  $: D_recov_hospital  = 21

  $: state = location.protocol + '//' + location.host + location.pathname + "?" + queryString.stringify({
   "Country":country,
   "logN":logN.toFixed(2),
   "I0":I0.toFixed(2),
   "R0":R0.toFixed(2),
   "InterventionTime":InterventionTime.toFixed(2),
   "InterventionAmt":InterventionAmt.toFixed(2),
   "N_ICU":N_ICU.toFixed(2),
   "P_symptoms": P_symptoms.toFixed(2),
   "P_hospital": P_hospital.toFixed(2),
   "P_ICU": P_ICU.toFixed(2),
   "P_die_in_ICU": P_die_in_ICU.toFixed(2),
   "P_die_in_hospital": P_die_in_hospital.toFixed(2),
   "D_to_asym": D_to_asym.toFixed(2),
   "D_to_symptoms": D_to_symptoms.toFixed(2),
   "D_to_hospital": D_to_hospital.toFixed(2),
   "D_to_ICU": D_to_ICU.toFixed(2),
   "D_to_death": D_to_death.toFixed(2),
   "D_recov_asym": D_recov_asym.toFixed(2),
   "D_recov_symptoms": D_recov_symptoms.toFixed(2),
   "D_recov_ICU": D_recov_ICU.toFixed(2),
   "D_recov_hospital": D_recov_hospital.toFixed(2),
  })


  function get_solution(dt, N, I0, R0, P_symptoms, P_hospital, P_ICU, P_die_in_ICU, P_die_in_hospital, D_to_asym, D_to_symptoms, D_to_hospital, D_to_ICU, D_to_death, D_recov_asym, D_recov_symptoms, D_recov_ICU, D_recov_hospital, N_ICU, InterventionTime, InterventionAmt) {

    var interpolation_steps = 40
    var steps = 110*interpolation_steps
    var dt = dt/interpolation_steps
    var sample_step = interpolation_steps

    var method = Integrators["RK4"]
    /*function f(t, x){

      // SEIR ODE
      if (t > InterventionTime && t < InterventionTime + duration){
        var beta = (InterventionAmt)*R0/(D_infectious)
      } else if (t > InterventionTime + duration) {
        var beta = 0.5*R0/(D_infectious)        
      } else {
        var beta = R0/(D_infectious)
      }
      var a     = 1/D_incubation
      var gamma = 1/D_infectious
      
      var S        = x[0] // Susectable
      var E        = x[1] // Exposed
      var I        = x[2] // Infectious 
      var Mild     = x[3] // Recovering (Mild)     
      var Severe   = x[4] // Recovering (Severe at home)
      var Severe_H = x[5] // Recovering (Severe in hospital)
      var Fatal    = x[6] // Recovering (Fatal)
      var R_Mild   = x[7] // Recovered
      var R_Severe = x[8] // Recovered
      var R_Fatal  = x[9] // Dead

      var p_severe = P_SEVERE
      var p_fatal  = CFR
      var p_mild   = 1 - P_SEVERE - CFR

      var dS        = -beta*I*S
      var dE        =  beta*I*S - a*E
      var dI        =  a*E - gamma*I
      var dMild     =  p_mild*gamma*I   - (1/D_recovery_mild)*Mild
      var dSevere   =  p_severe*gamma*I - (1/D_hospital_lag)*Severe
      var dSevere_H =  (1/D_hospital_lag)*Severe - (1/D_recovery_severe)*Severe_H
      var dFatal    =  p_fatal*gamma*I  - (1/D_death)*Fatal
      var dR_Mild   =  (1/D_recovery_mild)*Mild
      var dR_Severe =  (1/D_recovery_severe)*Severe_H
      var dR_Fatal  =  (1/D_death)*Fatal

      //      0   1   2   3      4        5          6       7        8          9
      return [dS, dE, dI, dMild, dSevere, dSevere_H, dFatal, dR_Mild, dR_Severe, dR_Fatal]
    }*/
    function f(t, x){

      // SEIR ODE
      if (t > InterventionTime){
        var beta = (InterventionAmt)*R0/(D_to_symptoms+D_recov_symptoms)
      } else {
        var beta = R0/(D_to_symptoms+D_recov_symptoms)
      }
      if (true){  // correction if days are median instead of average
		  var crct = Math.log(2)
	  } else {
		  var crct = 1
	  }
	  var eps = 1e-1
      
      var r_ea  = crct/D_to_asym
      
      var D_recov_asym = D_to_symptoms
      
      var r_ai  = crct/D_to_symptoms
      var r_aw  = crct/D_to_symptoms * (1-P_symptoms) / P_symptoms
      var r_wr  = crct/Math.max(eps, D_recov_asym - D_to_symptoms * P_symptoms / (1-P_symptoms))
      var r_ad  = 0
      
      var r_ih  = crct/D_to_hospital
      var r_ix  = crct/D_to_hospital * (1-P_hospital)/P_hospital
      var r_xr  = crct/Math.max(eps, D_recov_symptoms - D_to_hospital * P_hospital / (1-P_hospital))
      var r_id  = 0
      
      var r_hj  = crct/D_to_ICU
      var r_hy  = crct/D_to_ICU * (1-P_ICU)/P_ICU
      var r_yr  = crct/Math.max(eps, D_recov_hospital - D_to_ICU * P_ICU / (1-P_ICU))
      var r_hd  = 0
      
      var r_jd  = crct/D_to_death
      var r_jz  = crct/D_to_death * (1-P_die_in_ICU)/P_die_in_ICU
      var r_zr  = crct/Math.max(eps, D_recov_ICU - D_to_death * P_die_in_ICU / (1-P_die_in_ICU))
      
      var S        = x[0] // Susceptable
      var E        = x[1] // Exposed
      var A        = x[2] // Infectious 
      var W        = x[3] // Infectious 
      var I        = x[4] // Infectious 
      var X        = x[5] // Infectious 
      var H        = x[6] // Hospital
      var Y        = x[7] // Hospital
      var J        = x[8] // ICU
      var Z        = x[9] // ICU
      var R_A      = x[10] // Recovering (Mild)     
      var R_I      = x[11] // Recovering (Mild)     
      var R_H      = x[12] // Recovering (Hospital)
      var R_J      = x[13] // Recovering (Hospital)
      var D_A      = x[14] // Dead (Mild)
      var D_I      = x[15] // Dead (Mild)
      var D_H      = x[16] // Dead (Hospital)
      var D_J      = x[17] // Dead (ICU)


      var dS       = -beta*(A+I+H+J+W+X+Y+Z)*S
      var dE       = beta*(A+I+H+J+W+X+Y+Z)*S - r_ea * E
      var dA       = r_ea * E - (r_ai+r_aw+r_ad)*A
      var dW       = r_aw * A - r_wr*W
      var dI       = r_ai * A - (r_ih+r_ix+r_id)*I
      var dX       = r_ix * I - r_xr*X
      var dH       = r_ih * I - (r_hj+r_hy+r_hd)*H
      var dY       = r_hy * H - r_yr*Y
      var dJ       = r_hj * H - (r_jd+r_jz)*J
      var dZ       = r_jz * J - r_zr*Z
      
      var dR_A     = r_wr*W
      var dR_I     = r_xr*X
      var dR_H     = r_yr*Y
      var dR_J     = r_zr*Z
      
      var dD_A     = r_ad*A
      var dD_I     = r_id*I
      var dD_H     = r_hd*H
      var dD_J     = r_jd*J
      
      if (N*(J+Z) >= N_ICU*N/1e5) {
		if (dJ+dZ > 0){
			dD_H = dJ+dZ
			dJ = 0
			dZ = 0
		}	
	  }
	  
      //      0   1   2   3   4   5   6   7   8   9   10    11    12    13    14    15    16    17
      return [dS, dE, dA, dW, dI, dX, dH, dY, dJ, dZ, dR_A, dR_I, dR_H, dR_J, dD_A, dD_I, dD_H, dD_J]
    }

    var v = [1 - I0/N, 0, I0/N*(1-P_symptoms), 0, I0/N*P_symptoms, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    var t = 0

    var P  = []
    var TI = []
    var Iters = []
    var Plotting = []
    var dPlotting = []
    while (steps--) { 
      if ((steps+1) % (sample_step) == 0) {
            //    Dead                        ICU            Hospital       Recovered                    Sick           Infectious
        P.push([ N*(v[14]+v[15]+v[16]+v[17]), N*(v[8]+v[9]), N*(v[6]+v[7]), N*(v[10]+v[11]+v[12]+v[13]), N*(v[4]+v[5]), N*(v[2] + v[3]) ])
        Iters.push(v)
        TI.push(N*(1-v[0]))
        //console.log((v[0] + v[1] + v[2] + v[3] + v[4] + v[5] + v[6] + v[7] + v[8] + v[9] + v[10] + v[11] + v[12] + v[13] + v[14] + v[15] + v[16]+ v[17]))
        //console.log(v[0] , v[1] , v[2] , v[3] , v[4] , v[5] , v[6] , v[7] , v[8] , v[9] , v[10] , v[11] , v[12] , v[13] , v[14] , v[15] , v[16]+ v[17])
      }
      v =integrate(method,f,v,t,dt); 
      t+=dt
    }
    return {"P": P, 
            "deaths": N*(v[14]+v[15]+v[16]+v[17]), 
            "total": 1-v[0],
            "total_infected": TI,
            "Iters":Iters,
            "Plotting":Plotting,
            "dIters": f}
  }

  function max(P, checked) {
    return P.reduce((max, b) => Math.max(max, sum(b, checked) ), sum(P[0], checked) )
  }

  $: Sol            = get_solution(dt, N, I0, R0, P_symptoms, P_hospital, P_ICU, P_die_in_ICU, P_die_in_hospital, D_to_asym, D_to_symptoms, D_to_hospital, D_to_ICU, D_to_death, D_recov_asym, D_recov_symptoms, D_recov_ICU, D_recov_hospital, N_ICU, InterventionTime, InterventionAmt)
  $: P              = Sol["P"].slice(0,100)
  $: timestep       = dt
  $: tmax           = dt*100
  $: deaths         = Sol["deaths"]
  $: total          = Sol["total"]
  $: total_infected = Sol["total_infected"].slice(0,100)
  $: Iters          = Sol["Iters"]
  $: dIters         = Sol["dIters"]
  $: Plotting       = Sol["Plotting"]
  $: Pmax           = max(P, checked)
  $: lock           = false

  var colors = [ "#386cb0", "#aa0086", "#8da0cb", "#4daf4a", "#f0027f", "#fdc086"]

  var Plock = 1

  var drag_y = function (){
    var dragstarty = 0
    var Pmaxstart = 0

    var dragstarted = function (d) {
      dragstarty = event.y  
      Pmaxstart  = Pmax
    }

    var dragged = function (d) {
      Pmax = Math.max( (Pmaxstart*(1 + (event.y - dragstarty)/500)), 10)
    }

    return drag().on("drag", dragged).on("start", dragstarted)
  }

  var drag_x = function (){
    var dragstartx = 0
    var dtstart = 0
    var Pmaxstart = 0
    var dragstarted = function (d) {
      dragstartx = event.x
      dtstart  = dt
      Plock = Pmax
      lock = true
    }
    var dragged = function (d) {
      dt = dtstart - 0.0015*(event.x - dragstartx)
    }
    var dragend = function (d) {
      lock = false
    }
    return drag().on("drag", dragged).on("start", dragstarted).on("end", dragend)
  }

  var drag_intervention = function (){
    var dragstarty = 0
    var InterventionTimeStart = 0

    var dragstarted = function (d) {
      dragstarty = event.x  
      InterventionTimeStart = InterventionTime
      Plock = Pmax
      lock = true
    }

    var dragged = function (d) {
      // InterventionTime = Math.max( (*(1 + (event.x - dragstarty)/500)), 10)
      // console.log(event.x)
      InterventionTime = Math.min(tmax-1, Math.max(0, InterventionTimeStart + xScaleTimeInv(event.x - dragstarty)))
    }

    var dragend = function (d) {
      lock = false
    }

    return drag().on("drag", dragged).on("start", dragstarted).on("end", dragend)
  }


  var drag_intervention_end = function (){
    var dragstarty = 0
    var durationStart = 0

    var dragstarted = function (d) {
      dragstarty = event.x  
      durationStart = duration
      Plock = Pmax
      lock = true
    }

    var dragged = function (d) {
      // InterventionTime = Math.max( (*(1 + (event.x - dragstarty)/500)), 10)
      // console.log(event.x)
      duration = Math.min(tmax-1, Math.max(0, durationStart + xScaleTimeInv(event.x - dragstarty)))
    }

    var dragend = function (d) {
      lock = false
    }

    return drag().on("drag", dragged).on("start", dragstarted).on("end", dragend)
  }


  $: parsed = "";
  onMount(async () => {
    var drag_callback_y = drag_y()
    drag_callback_y(selectAll("#yAxisDrag"))
    var drag_callback_x = drag_x()
    drag_callback_x(selectAll("#xAxisDrag"))
    var drag_callback_intervention = drag_intervention()
    // drag_callback_intervention(selectAll("#interventionDrag"))
    drag_callback_intervention(selectAll("#dottedline"))
    // var drag_callback_intervention_end = drag_intervention_end()
    // drag_callback_intervention_end(selectAll("#dottedline2"))

    if (typeof window !== 'undefined') {
      parsed = queryString.parse(window.location.search)
      if (!(parsed.logN === undefined)) {logN = parseFloat(parsed.logN)}
      if (!(parsed.I0 === undefined)) {I0 = parseFloat(parsed.I0)}
      if (!(parsed.R0 === undefined)) {R0 = parseFloat(parsed.R0)}
      if (!(parsed.InterventionTime === undefined)) {InterventionTime = parseFloat(parsed.InterventionTime)}
      if (!(parsed.InterventionAmt === undefined)) {InterventionAmt = parseFloat(parsed.InterventionAmt)}
      if (!(parsed.N_ICU === undefined)) {N_ICU = parseFloat(parsed.N_ICU)}
      if (!(parsed.P_symptoms === undefined)) {P_symptoms = parseFloat(parsed.P_symptoms)}
      if (!(parsed.P_hospital === undefined)) {P_hospital = parseFloat(parsed.P_hospital)}
      if (!(parsed.P_ICU === undefined)) {P_ICU = parseFloat(parsed.P_ICU)}
      if (!(parsed.P_die_in_ICU === undefined)) {P_die_in_ICU = parseFloat(parsed.P_die_in_ICU)}
      if (!(parsed.P_die_in_hospital === undefined)) {P_die_in_hospital = parseFloat(parsed.P_die_in_hospital)}
      if (!(parsed.D_to_asym === undefined)) {D_to_asym = parseFloat(parsed.D_to_asym)}
      if (!(parsed.D_to_symptoms === undefined)) {D_to_symptoms = parseFloat(parsed.D_to_symptoms)}
      if (!(parsed.D_to_hospital === undefined)) {D_to_hospital = parseFloat(parsed.D_to_hospital)}
      if (!(parsed.D_to_ICU === undefined)) {D_to_ICU = parseFloat(parsed.D_to_ICU)}
      if (!(parsed.D_to_death === undefined)) {D_to_death = parseFloat(parsed.D_to_death)}
      if (!(parsed.D_recov_asym === undefined)) {D_recov_asym = parseFloat(parsed.D_recov_asym)}
      if (!(parsed.D_recov_symptoms === undefined)) {D_recov_symptoms = parseFloat(parsed.D_recov_symptoms)}
      if (!(parsed.D_recov_ICU === undefined)) {D_recov_ICU = parseFloat(parsed.D_recov_ICU)}
      if (!(parsed.D_recov_hospital === undefined)) {D_recov_hospital = parseFloat(parsed.D_recov_hospital)}
      if (!(parsed.country === undefined)) {country = parsed.country}
    }
  });

  function lock_yaxis(){
    Plock = Pmax
    lock  = true
  }

  function unlock_yaxis(){
    lock = false
  }

  const padding = { top: 20, right: 0, bottom: 20, left: 25 };

  let width  = 750;
  let height = 400;

  $: xScaleTime = scaleLinear()
    .domain([0, tmax])
    .range([padding.left, width - padding.right]);

  $: xScaleTimeInv = scaleLinear()
    .domain([0, width])
    .range([0, tmax]);

  $: indexToTime = scaleLinear()
    .domain([0, P.length])
    .range([0, tmax])

  window.addEventListener('mouseup', unlock_yaxis);

  $: checked = [true, true, true, false, true, true]
  $: active  = 0
  $: active_ = active >= 0 ? active : Iters.length - 1

  var Tinc_s = "\\color{#CCC}{T^{-1}_{\\text{inc}}} "
  var Tinf_s = "\\color{#CCC}{T^{-1}_{\\text{inf}}}"
  var Rt_s   = "\\color{#CCC}{\\frac{\\mathcal{R}_{t}}{T_{\\text{inf}}}} "
  $: ode_eqn = katex.renderToString("\\frac{d S}{d t}=-" +Rt_s +"\\cdot IS,\\qquad \\frac{d E}{d t}=" +Rt_s +"\\cdot IS- " + Tinc_s + " E,\\qquad \\frac{d I}{d t}=" + Tinc_s + "E-" + Tinf_s+ "I, \\qquad \\frac{d R}{d t}=" + Tinf_s+ "I", {
    throwOnError: false,
    displayMode: true,
    colorIsTextColor: true
  });

  function math_inline(str) {
    return katex.renderToString(str, {
    throwOnError: false,
    displayMode: false,
    colorIsTextColor: true
    });
  }

  function math_display(str) {
    return katex.renderToString(str, {
    throwOnError: false,
    displayMode: true,
    colorIsTextColor: true
    });
  }
  
  $: p_num_ind = 40

  $: get_d = function(i){
    return dIters(indexToTime(i), Iters[i])
  }

  function get_milestones(P){

    function argmax(x, index) {
      return x.map((x, i) => [x[index], i]).reduce((r, a) => (a[0] > r[0] ? a : r))[1];
    }

	//    Dead                        ICU            Hospital       Recovered                    Sick           Infectious
    var milestones = []
    for (var i = 0; i < P.length; i++) {
      if (P[i][0] >= 0.5) {
        milestones.push([i*dt, "First death"])
        break
      }
    }

    var i = argmax(P, 2)
    milestones.push([i*dt, "</br>Peak hospital: " + format(",")(Math.round(P[i][2]))])
    var j = argmax(P, 1)
    milestones.push([j*dt, "Peak ICU: " + format(",")(Math.round(P[j][1]))])

    return milestones
  }

  $: milestones = get_milestones(P)
  $: log = true

</script>

<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.11.1/dist/katex.css" integrity="sha384-bsHo4/LA+lkZv61JspMDQB9QP1TtO4IgOf2yYS+J6VdAYLVyx1c3XKcsHh0Vy8Ws" crossorigin="anonymous">

<style>
  .small { font: italic 6px Source Code Pro; }
  @import url('https://fonts.googleapis.com/css?family=Source+Code+Pro&display=swap');


  h2 {
    margin: auto;
    width: 950px;
    font-size: 40px;
    padding-top: 20px;
    padding-bottom: 20px;
    font-weight: 300;
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    padding-bottom: 30px
  }

  .center {
    margin: auto;
    width: 950px;
    padding-bottom: 20px;
    font-weight: 300;
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    color:#666;
    font-size: 16.5px;
    text-align: justify;
    line-height: 24px
  }

  .ack {
    margin: auto;
    width: 950px;
    padding-bottom: 20px;
    font-weight: 300;
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    color:#333;
    font-size: 13px;
  }

  .row {
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    margin: auto;
    display: flex;
    width: 948px;
    font-size: 13px;
  }

  .caption {
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    font-size: 13px;    
  }

  .column {
    flex: 158px;
    padding: 0px 5px 5px 0px;
    margin: 0px 5px 5px 5px;
    /*border-top: 2px solid #999*/
  }

  .minorTitle {
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    margin: auto;
    display: flex;
    width: 950px;
    font-size: 17px;
    color: #666;
  }

  .minorTitleColumn{
    flex: 60px;
    padding: 3px;
    border-bottom: 2px solid #999;
  }


  .paneltext{
    position:relative;
    height:130px;
  }

  .paneltitle{
    color:#777; 
    line-height: 17px; 
    padding-bottom: 4px;
    font-weight: 700;
    font-family: nyt-franklin,helvetica,arial,sans-serif;
  }

  .paneldesc{
    color:#888; 
    text-align: left;
    font-weight: 300;
  }

  .slidertext{
    color:#555; 
    line-height: 7px; 
    padding-bottom: 0px; 
    padding-top: 7px;
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    font-family: 'Source Code Pro', monospace;
    font-size: 10px;
    text-align: right;
    /*font-weight: bold*/
  }
    
  .range {
    width: 100%;
  }

  .chart {
    width: 100%;
    margin: 0 auto;
    padding-top:0px;
    padding-bottom:10px;
  }

  .legend {
    color: #888;
    font-family: Helvetica, Arial;
    font-size: .725em;
    font-weight: 200;
    height: 100px;
    left: 20px;
    top: 4px;
    position: absolute;
  }

  .legendtitle {
    color:#777; 
    font-size:13px;
    padding-bottom: 6px;
    font-weight: 600;
    font-family: nyt-franklin,helvetica,arial,sans-serif;
  }


  .legendtext{
    color:#888; 
    font-size:13px;
    padding-bottom: 5px;
    font-weight: 300;
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    line-height: 14px;
  }

  .legendtextnum{
    color:#888; 
    font-size:13px;
    padding-bottom: 5px;
    font-weight: 300;
    line-height: 12px;
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    left: -3px;
    position: relative;
  }

  .tick {
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    font-size: .725em;
    font-weight: 200;
    font-size: 13px
  }

  td { 
    text-align: left;
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    border-bottom: 1px solid #DDD;
    border-collapse: collapse;
    padding: 3px;
    /*font-size: 14px;*/
  }

  tr {
    border-collapse: collapse;
    border-spacing: 15px;
  }

  .eqn {
    font-family: nyt-franklin,helvetica,arial,sans-serif;
    margin: auto;
    display: flex;
    flex-flow: row wrap;
    width: 950px;
    column-count: 4;
    font-weight: 300;
    color:#666;
    font-size: 16.5px;
  }

  th { font-weight: 500; text-align: left; padding-bottom: 5px; vertical-align: text-top;     border-bottom: 1px solid #DDD; }

  a:link { color: grey; }
  a:visited { color: grey; }

</style>

<h2>Epidemic Calculator (adapted for {country})</h2>

<div class="chart" style="display: flex; max-width: 1120px">

  <div style="flex: 0 0 290px; width:270px;">
    <div style="position:relative; top:48px; right:-115px">
      <div class="legendtext" style="position:absolute; left:-16px; top:-34px; width:50px; height: 100px; font-size: 13px; line-height:16px; font-weight: normal; text-align: center"><b>Day</b><br> {Math.round(indexToTime(active_))}</div>

      <!-- Susceptible -->
      <div style="position:absolute; left:0px; top:0px; width: 180px; height: 100px">

        <span style="pointer-events: none"><Checkbox color="#CCC"/></span>
        <Arrow height="41"/>

        <div class="legend" style="position:absolute;">
          <div class="legendtitle">Susceptible</div>
          <div style="padding-top: 5px; padding-bottom: 1px">
          <div class="legendtextnum"><span style="font-size:12px; padding-right:3px; color:#CCC">∑</span> <i>{formatNumber(Math.round(N*Iters[active_][0]))} 
                                  ({ (100*Iters[active_][0]).toFixed(2) }%)</i></div>
          <div class="legendtextnum"><span style="font-size:12px; padding-right:2px; color:#CCC">Δ</span> <i>{formatNumber(Math.round(N*get_d(active_)[0]))} / day</i>
                                 </div>
          </div>
        </div>
          <div class="legendtext" style="text-align: right; width:105px; left:-111px; top: 4px; position:relative;">Population not immune to disease</div>

      </div>

      <!-- Infectious -->
      <div style="position:absolute; left:0px; top:{legendheight*1}px; width: 180px; height: 100px">

        <Checkbox color="{colors[5]}" bind:checked={checked[5]}/>      
        <Arrow height="41"/>

        <div class="legend" style="position:absolute;">
          <div class="legendtitle">Infectious</div>

          <div style="padding-top: 5px; padding-bottom: 1px">
          <div class="legendtextnum"><span style="font-size:12px; padding-right:3px; color:#CCC">∑</span> <i>{formatNumber(Math.round(N*(Iters[active_][2]+Iters[active_][3])))} 
                                  ({ (100*(Iters[active_][2]+Iters[active_][3])).toFixed(2) }%)</div>
          <div class="legendtextnum"><span style="font-size:12px; padding-right:2px; color:#CCC">Δ</span> <i>{formatNumber(Math.round(N*(get_d(active_)[2]+get_d(active_)[3]))) } / day</i>
                                 </div>
          </div>
        </div>
        <div class="legendtext" style="text-align: right; width:105px; left:-111px; top: 4px; position:relative;">People feeling fine but infecting others</div>

      </div>

      <!-- Sick -->
      <div style="position:absolute; left:0px; top:{legendheight*2}px; width: 180px; height: 100px">

        <Checkbox color="{colors[4]}" bind:checked={checked[4]}/>
        <Arrow height="41"/>   

        <div class="legend" style="position:absolute;">
          <div class="legendtitle">Sick</div>
          <div style="padding-top: 5px; padding-bottom: 1px">
          <div class="legendtextnum"><span style="font-size:12px; padding-right:3px; color:#CCC">∑</span> <i>{formatNumber(Math.round(N*(Iters[active_][4]+Iters[active_][5])))}
                                  ({ (100*(Iters[active_][4]+Iters[active_][5])).toFixed(2) }%)</div>
          <div class="legendtextnum"><span style="font-size:12px; padding-right:2px; color:#CCC">Δ</span> <i>{formatNumber(Math.round(N*(get_d(active_)[4]+get_d(active_)[5]))) } / day</i>
                                 </div>
          </div>
        </div>
        <div class="legendtext" style="text-align: right; width:105px; left:-111px; top: 4px; position:relative;">People feeling <i>sick</i></div>


      </div>

      <!-- Removed -->
      <div style="position:absolute; left:0px; top:{legendheight*3}px; width: 180px; height: 100px">

        <Checkbox color="grey" callback={(s) => {checked[0] = s; checked[1] = s; checked[2] = s; checked[3] = s} }/>
        <Arrow height="56" arrowhead="" dasharray="3 2"/>

        <div class="legend" style="position:absolute;">
          <div class="legendtitle">Removed</div>
          <div style="padding-top: 10px; padding-bottom: 1px">
          <div class="legendtextnum"><span style="font-size:12px; padding-right:3px; color:#CCC">∑</span> <i>{formatNumber(Math.round(N* (Iters[active_][10]+Iters[active_][11]+Iters[active_][1]+Iters[active_][12]+Iters[active_][13]+Iters[active_][14]+Iters[active_][15]+Iters[active_][16]+Iters[active_][17]) ))} 
                                  ({ ((100*(Iters[active_][10]+Iters[active_][11]+Iters[active_][12]+Iters[active_][13]+Iters[active_][14]+Iters[active_][15]+Iters[active_][16]+Iters[active_][17]))).toFixed(2) }%)</div>
          <div class="legendtextnum"><span style="font-size:12px; padding-right:2px; color:#CCC">Δ</span> <i>{formatNumber(Math.round(N*(get_d(active_)[10]+get_d(active_)[11]+get_d(active_)[12]+get_d(active_)[13]+get_d(active_)[14]+get_d(active_)[15]+get_d(active_)[16]+get_d(active_)[17]) )) } / day</i>
                                 </div>
          </div>
        </div>
        <div class="legendtext" style="text-align: right; width:105px; left:-111px; top: 4x; position:relative;">Population no longer infectious due to isolation or immunity</div>

      </div>

      <!-- Recovered -->
      <div style="position:absolute; left:0px; top:{legendheight*4+14-3}px; width: 180px; height: 100px">
        <Checkbox color="{colors[3]}" bind:checked={checked[3]}/>
        <Arrow height="23" arrowhead="" dasharray="3 2"/>
        <div class="legend" style="position:absolute;">
          <div class="legendtitle">Recovered</div>

          <div style="padding-top: 3px; padding-bottom: 1px">
          <div class="legendtextnum"><span style="font-size:12px; padding-right:3px; color:#CCC">∑</span> <i>{formatNumber(Math.round(N*(Iters[active_][10]+Iters[active_][11]+Iters[active_][12]+Iters[active_][13]) ))} 
                                  ({ (100*(Iters[active_][10]+Iters[active_][11]+Iters[active_][12]+Iters[active_][13])).toFixed(2) }%)</div>
          </div>
        </div>
        <div class="legendtext" style="text-align: right; width:105px; left:-111px; top: 8px; position:relative;">Full recoveries</div>

      </div>

      <!-- Hospitalized -->
      <div style="position:absolute; left:0px; top:{legendheight*4+57}px; width: 180px; height: 100px">
        <Arrow height="43" arrowhead="" dasharray="3 2"/>
        <Checkbox color="{colors[2]}" bind:checked={checked[2]}/>
        <div class="legend" style="position:absolute;">
          <div class="legendtitle">Hospitalized</div>
          <div style="padding-top: 3px; padding-bottom: 1px">
          <div class="legendtextnum"><span style="font-size:12px; padding-right:3px; color:#CCC">∑</span> <i>{formatNumber(Math.round(N*(Iters[active_][6]+Iters[active_][7]) ))} 
                                  ({ (100*(Iters[active_][6]+Iters[active_][7])).toFixed(2) }%)</div>
          </div>
          <div class="legendtextnum"><span style="font-size:12px; padding-right:2px; color:#CCC">Δ</span> <i>{formatNumber(Math.round(N*(get_d(active_)[6]+get_d(active_)[7]))) } / day</i>
                                 </div>
        </div>
        <div class="legendtext" style="text-align: right; width:105px; left:-111px; top: 10px; position:relative;">Active hospitalizations</div>

      </div>

       <!-- ICU -->
      <div style="position:absolute; left:0px; top:{legendheight*4+120+2}px; width: 180px; height: 100px">
        <Arrow height="43" arrowhead="" dasharray="3 2"/>
        <Checkbox color="{colors[1]}" bind:checked={checked[1]}/>
        <div class="legend" style="position:absolute;">
          <div class="legendtitle">ICU</div>
          <div style="padding-top: 3px; padding-bottom: 1px">
          <div class="legendtextnum"><span style="font-size:12px; padding-right:3px; color:#CCC">∑</span> <i>{formatNumber(Math.round(N*(Iters[active_][8]+Iters[active_][9]) ))} 
                                  ({ (100*(Iters[active_][8]+Iters[active_][9])).toFixed(2) }%)</div>
          </div>
          <div class="legendtextnum"><span style="font-size:12px; padding-right:2px; color:#CCC">Δ</span> <i>{formatNumber(Math.round(N*(get_d(active_)[8]+get_d(active_)[9]))) } / day</i>
                                 </div>
        </div>
        <div class="legendtext" style="text-align: right; width:105px; left:-111px; top: 10px; position:relative;">People in the Intensive Care Unit</div>

      </div>

      <div style="position:absolute; left:0px; top:{legendheight*4+120+2+63+2}px; width: 180px; height: 100px">
        <Arrow height="40" arrowhead="" dasharray="3 2"/>

        <Checkbox color="{colors[0]}" bind:checked={checked[0]}/>

        <div class="legend" style="position:absolute;">
          <div class="legendtitle">Fatalities</div>
          <div style="padding-top: 3px; padding-bottom: 1px">          
          <div class="legendtextnum"><span style="font-size:12px; padding-right:3px; color:#CCC">∑</span> <i>{formatNumber(Math.round(N*(Iters[active_][14]+Iters[active_][15]+Iters[active_][16]+Iters[active_][17])))} 
                                  ({ (100*(Iters[active_][14]+Iters[active_][15]+Iters[active_][16]+Iters[active_][17])).toFixed(2) }%)</div>
          <div class="legendtextnum"><span style="font-size:12px; padding-right:2px; color:#CCC">Δ</span> <i>{formatNumber(Math.round(N*(get_d(active_)[14]+get_d(active_)[15]+get_d(active_)[16]+get_d(active_)[17]))) } / day</i>
                                 </div>
          </div>
        </div>
        <div class="legendtext" style="text-align: right; width:105px; left:-111px; top: 10px; position:relative;">Deaths</div>
      </div>
    </div>
  </div>

  <div style="flex: 0 0 890px; width:890px; height: {height+128+60}px; position:relative;">

    <div style="position:relative; top:60px; left: 10px">
      <Chart bind:checked={checked}
             bind:active={active}
             y = {P} 
             xmax = {Xmax} 
             total_infected = {total_infected} 
             deaths = {deaths} 
             total = {total} 
             timestep={timestep}
             tmax={tmax}
             N={N}
             ymax={lock ? Plock: Pmax}
             InterventionTime={InterventionTime}
             colors={colors}
             log={!log}/>
      </div>

      <div id="xAxisDrag"
           style="pointer-events: all;
                  position: absolute;
                  top:{height+80}px;
                  left:{0}px;
                  width:{780}px;
                  background-color:#222;
                  opacity: 0;
                  height:25px;
                  cursor:col-resize">
      </div>

      <div id="yAxisDrag"
           style="pointer-events: all;
                  position: absolute;
                  top:{55}px;
                  left:{0}px;
                  width:{20}px;
                  background-color:#222;
                  opacity: 0;
                  height:425px;
                  cursor:row-resize">
      </div>

      <!-- Intervention Line -->
      <div style="position: absolute; width:{width+15}px; height: {height}px; position: absolute; top:100px; left:10px; pointer-events: none">
        <div id="dottedline"  style="pointer-events: all;
                    position: absolute;
                    top:-38px;
                    left:{xScaleTime(InterventionTime)}px;
                    visibility: {(xScaleTime(InterventionTime) < (width - padding.right)) ? 'visible':'hidden'};
                    width:2px;
                    background-color:#FFF;
                    border-right: 1px dashed black;
                    pointer-events: all;
                    cursor:col-resize;
                    height:{height+19}px">

        <div style="position:absolute; opacity: 0.5; top:-5px; left:10px; width: 120px">
        <span style="font-size: 13px">{@html math_inline("\\mathcal{R}_t=" + (R0*InterventionAmt).toFixed(2) )}</span> ⟶ 
        </div>

        {#if xScaleTime(InterventionTime) >= 100}
          <div style="position:absolute; opacity: 0.5; top:-2px; left:-97px; width: 120px">
          <span style="font-size: 13px">⟵ {@html math_inline("\\mathcal{R}_0=" + (R0).toFixed(2) )}</span>
          </div>      
        {/if}

        <div id="interventionDrag" class="legendtext" style="flex: 0 0 160px; width:120px; position:relative;  top:-70px; height: 60px; padding-right: 15px; left: -125px; pointer-events: all;cursor:col-resize;" >
          <div class="paneltitle" style="top:9px; position: relative; text-align: right">Intervention on day {format("d")(InterventionTime)}</div>
          <span></span><div style="top:9px; position: relative; text-align: right">
          (drag me)</div>
          <div style="top:43px; left:40px; position: absolute; text-align: right; width: 20px; height:20px; opacity: 0.3">
            <svg width="20" height="20">
              <g transform="rotate(90)">
                <g transform="translate(0,-20)">
                  <path d="M2 11h16v2H2zm0-4h16v2H2zm8 11l3-3H7l3 3zm0-16L7 5h6l-3-3z"/>
                 </g>  
              </g>
            </svg>
          </div>
        </div>


        <div style="width:150px; position:relative; top:-85px; height: 80px; padding-right: 15px; left: 0px; ;cursor:col-resize; background-color: white; position:absolute" >

        </div>


        </div>
      </div>

      <!-- Intervention Line slider -->
      <div style="position: absolute; width:{width+15}px; height: {height}px; position: absolute; top:120px; left:10px; pointer-events: none">
        <div style="
            position: absolute;
            top:-38px;
            left:{xScaleTime(InterventionTime)}px;
            visibility: {(xScaleTime(InterventionTime) < (width - padding.right)) ? 'visible':'hidden'};
            width:2px;
            background-color:#FFF;
            border-right: 1px dashed black;
            cursor:col-resize;
            height:{height}px">
            <div style="flex: 0 0 160px; width:200px; position:relative; top:-125px; left: 1px" >
              <div class="caption" style="pointer-events: none; position: absolute; left:0; top:40px; width:150px; border-left: 2px solid #777; padding: 5px 7px 7px 7px; ">      
              <div class="paneltext"  style="height:20px; text-align: right">
              <div class="paneldesc">to decrease transmission by<br></div>
              </div>
              <div style="pointer-events: all">
              <div class="slidertext" on:mousedown={lock_yaxis}>{(100*(1-InterventionAmt)).toFixed(2)}%</div>
              <input class="range" type=range bind:value={OMInterventionAmt} min=0 max=1 step=0.01 on:mousedown={lock_yaxis}>
              </div>
              </div>
            </div>
          </div>
      </div>

      <div style="pointer-events: none;
                  position: absolute;
                  top:{height+84}px;
                  left:{0}px;
                  width:{780}px;
                  opacity: 1.0;
                  height:25px;
                  cursor:col-resize">
            {#each milestones as milestone}
              <div style="position:absolute; left: {xScaleTime(milestone[0])+8}px; top: -30px;">
                <span style="opacity: 0.3"><Arrow height=30 arrowhead="#circle" dasharray = "2 1"/></span>
                  <div class="tick" style="position: relative; left: 0px; top: 35px; max-width: 130px; color: #BBB; background-color: white; padding-left: 4px; padding-right: 4px">{@html milestone[1]}</div>
              </div>
            {/each}
      </div>
    
    <div style="opacity:{xScaleTime(InterventionTime) >= 192? 1.0 : 0.2}">
      <div class="tick" style="color: #AAA; position:absolute; pointer-events:all; left:10px; top: 10px">
        <Checkbox color="#CCC" bind:checked={log}/><div style="position: relative; top: 4px; left:20px">linear scale</div>
      </div>
    </div>

   </div>

</div>



<div style="height:290px;">
  <div class="minorTitle">
    <div style="margin: 0px 0px 5px 4px; width:420px" class="minorTitleColumn">Transmission Dynamics</div>
    <div style="flex: 0 0 20; width:20px"></div>
    <div style="margin: 0px 4px 5px 0px" class="minorTitleColumn">Clinical Dynamics</div>
  </div>
  <div class = "row">

    <div class="column">
      <div class="paneltitle">Population Inputs</div>
      <div class="paneldesc" style="height:30px">Size of population<br></div>
      <div class="slidertext">{format(",")(Math.round(N))}</div>
      <input class="range" style="margin-bottom: 8px"type=range bind:value={logN} min={5} max=25 step=0.01>
      <div class="paneldesc" style="height:29px; border-top: 1px solid #EEE; padding-top: 10px">Number of initial infections<br></div>
      <div class="slidertext">{I0}</div>
      <input class="range" type=range bind:value={I0} min={1} max=1000 step=1>
    </div>

    <div class="column">
      <div class="paneltext">
      <div class="paneltitle">Basic Reproduction Number {@html math_inline("\\mathcal{R}_0")} </div>
      <div class="paneldesc">Measure of contagiousness: number of new infections per infected <br></div>
      </div>
      <div class="slidertext">{R0}</div>
      <input class="range" type=range bind:value={R0} min=0.01 max=10 step=0.01> 
      <div class="paneldesc" style="height:29px; border-top: 1px solid #EEE; padding-top: 10px">Days to being infectious<br></div>
      <div class="slidertext">{D_to_asym.toFixed(2)} Days</div>
      <input class="range" type=range bind:value={D_to_asym} min={0.1} max=30 step=0.1>      
    </div> 

    <div class="column">
      <div class="paneltitle">Incubation statistics</div>
      <div class="paneldesc" style="height:30px">Probability of staying asympromatic<br></div>
      <div class="slidertext">{((1-P_symptoms)*100).toFixed(2)} %</div>
      <input class="range" style="margin-bottom: 8px"type=range bind:value={P_symptoms} min=0 max=1 step=0.0001>
      <div class="paneldesc" style="height:29px; border-top: 1px solid #EEE; padding-top: 10px">Days to first symptoms<br></div>
      <div class="slidertext">{(D_to_asym+D_to_symptoms).toFixed(2)} Days</div>
      <input class="range" type=range bind:value={D_to_symptoms} min={0.0} max=30 step=0.1>
      <div class="paneldesc" style="height:29px; border-top: 1px solid #EEE; padding-top: 10px">Days to heal<br></div>
      <div class="slidertext">{D_recov_symptoms.toFixed(2)} Days</div>
      <input class="range" type=range bind:value={D_recov_symptoms} min={D_to_hospital * P_hospital / (1-P_hospital)} max=30 step=0.1>
    </div>

    <div class="column">
      <div class="paneltitle">Progression statistics</div>
      <div class="paneldesc" style="height:30px">Hospitalization probability<br></div>
      <div class="slidertext">{((P_hospital)*100).toFixed(2)} %</div>
      <input class="range" style="margin-bottom: 8px"type=range bind:value={P_hospital} min=0 max=1 step=0.0001>
      <div class="paneldesc" style="height:29px; border-top: 1px solid #EEE; padding-top: 10px">Days to get in hospital<br></div>
      <div class="slidertext">{D_to_hospital.toFixed(2)} Days</div>
      <input class="range" type=range bind:value={D_to_hospital} min={0.1} max=30 step=0.1>
      <div class="paneldesc" style="height:29px; border-top: 1px solid #EEE; padding-top: 10px">Days to recover from hospital<br></div>
      <div class="slidertext">{D_recov_hospital.toFixed(2)} Days</div>
      <input class="range" type=range bind:value={D_recov_hospital} min={D_to_ICU * P_ICU / (1-P_ICU)} max=30 step=0.1>
    </div>

    <div class="column">
      <div class="paneltitle">ICU statistics</div>
      <div class="paneldesc" style="height:30px">ICU probability from hospital<br></div>
      <div class="slidertext">{((P_ICU)*100).toFixed(2)} %</div>
      <input class="range" style="margin-bottom: 8px"type=range bind:value={P_ICU} min=0 max=1 step=0.0001>
      <div class="paneldesc" style="height:29px; border-top: 1px solid #EEE; padding-top: 10px">Days to get in ICU<br></div>
      <div class="slidertext">{D_to_ICU.toFixed(2)} Days</div>
      <input class="range" type=range bind:value={D_to_ICU} min={0.1} max=30 step=0.1>
      <div class="paneldesc" style="height:29px; border-top: 1px solid #EEE; padding-top: 10px">Days to recover from ICU<br></div>
      <div class="slidertext">{D_recov_ICU.toFixed(2)} Days</div>
      <input class="range" type=range bind:value={D_recov_ICU} min={D_to_death * P_die_in_ICU / (1-P_die_in_ICU)} max=30 step=0.1>
    </div>

    <div class="column">
      <div class="paneltitle">Mortality</div>
      <div class="paneldesc" style="height:30px">Probability to die in ICU<br></div>
      <div class="slidertext">{((P_die_in_ICU)*100).toFixed(2)} %</div>
      <input class="range" style="margin-bottom: 8px"type=range bind:value={P_die_in_ICU} min=0 max=1 step=0.0001>
      <div class="paneldesc" style="height:29px; border-top: 1px solid #EEE; padding-top: 10px">Days to die in ICU<br></div>
      <div class="slidertext">{D_to_death.toFixed(2)} Days</div>
      <input class="range" type=range bind:value={D_to_death} min={0.1} max=30 step=0.1>
      <div class="paneldesc" style="height:29px; border-top: 1px solid #EEE; padding-top: 10px">ICU beds per 100,000 inhabitants<br></div>
      <div class="slidertext">{N_ICU.toFixed(2)} beds</div>
      <input class="range" type=range bind:value={N_ICU} min={0.0} max=50 step=0.1>
    </div>

  </div>
</div>

<div style="position: relative; height: 12px"></div>

<p class = "center">
At the time of writing, the coronavirus disease of 2019 remains a global health crisis of grave and uncertain magnitude. To the non-expert (such as myself), contextualizing the numbers, forecasts and epidemiological parameters described in the media and literature can be challenging. I created this calculator as an attempt to address this gap in understanding.
</p>

<p class = "center">
A sampling of the estimates for epidemic parameters are presented below:
</p>

<div class="center">
<table style="width:100%; margin:auto; font-weight: 300; border-spacing: inherit">
  <tr>
    <th></th>
    <th>Location</th>
    <th>Reproduction Number<br> {@html math_inline("\\mathcal{R}_0")}</th>
    <th>Incubation Period<br> {@html math_inline("T_{\\text{inc}}")} (in days)</th>
    <th>Infectious Period<br> {@html math_inline("T_{\\text{inf}}")} (in days)</th>
  </tr>
  <tr>
    <td width="27%"><a href = "https://cmmid.github.io/topics/covid19/current-patterns-transmission/wuhan-early-dynamics.html">Kucharski et. al</a></td>
    <td>Wuhan </td>    
    <td>3.0 (1.5 — 4.5)</td>
    <td>5.2</td>
    <td>2.9</td>
  </tr>
  <tr>
    <td><a href = "https://www.nejm.org/doi/full/10.1056/NEJMoa2001316">Li, Leung and Leung</a></td>
    <td>Wuhan </td>    
    <td>2.2 (1.4 — 3.9)</td>
    <td>5.2 (4.1 — 7.0)</td>
    <td>2.3 (0.0 — 14.9)</td>
  </tr>
  <tr>
    <td><a href = "https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)30260-9/fulltext">Wu et. al</a></td>
    <td>Greater Wuhan </td>    
    <td>2.68 (2.47 — 2.86)</td>
    <td>6.1</td>
    <td>2.3</td>
  </tr>
  <tr>
    <td><a href = "https://www.who.int/news-room/detail/23-01-2020-statement-on-the-meeting-of-the-international-health-regulations-(2005)-emergency-committee-regarding-the-outbreak-of-novel-coronavirus-(2019-ncov)">WHO Initial Estimate</a></td>
    <td>Hubei </td>    
    <td>1.95 (1.4 — 2.5)</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td><a href = "https://www.who.int/docs/default-source/coronaviruse/who-china-joint-mission-on-covid-19-final-report.pdf">WHO-China Joint Mission </a></td>
    <td>Hubei </td>    
    <td>2.25 (2.0 — 2.5)</td>
    <td>5.5 (5.0 - 6.0)</td>
    <td></td>
  </tr>
  <tr>
    <td><a href = "https://www.biorxiv.org/content/10.1101/2020.01.25.919787v2">Liu et. al </a></td>
    <td>Guangdong</td>
    <td>4.5 (4.4 — 4.6)</td>
    <td>4.8 (2.2 — 7.4) </td>
    <td>2.9 (0 — 5.9)</td>
  </tr>
  <tr>
    <td><a href = "https://academic.oup.com/jtm/advance-article/doi/10.1093/jtm/taaa030/5766334">Rocklöv, Sjödin and Wilder-Smith</a></td>
    <td>Princess Diamond</td>
    <td>14.8</td>
    <td>5.0</td>
    <td>10.0</td>
  </tr>
  <tr>
    <td><a href = "https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.5.2000062">Backer, Klinkenberg, Wallinga</a></td>
    <td>Wuhan</td>
    <td></td>
    <td>6.5 (5.6 — 7.9)</td>
    <td></td>
  </tr>
  <tr>
    <td><a href = "https://www.medrxiv.org/content/10.1101/2020.01.23.20018549v2.article-info">Read et. al</a></td>
    <td>Wuhan</td>
    <td>3.11 (2.39 — 4.13)</td>
    <td></td>
    <td></td>
  </tr>
  <tr>
    <td><a href = "https://www.medrxiv.org/content/10.1101/2020.03.03.20028423v1">Bi et. al</a></td>
    <td>Shenzhen</td>
    <td></td>
    <td>4.8 (4.2 — 5.4)</td>
    <td>1.5 (0 — 3.4)</td>
    <td></td>
  </tr>

  <tr>
    <td><a href = "https://www.mdpi.com/2077-0383/9/2/462">Tang et. al</a></td>
    <td>China</td>
    <td>6.47 (5.71 — 7.23)</td>
    <td></td>
    <td></td>
  </tr>

</table>
</div>


<p class="center">
See [<a href="https://academic.oup.com/jtm/advance-article/doi/10.1093/jtm/taaa021/5735319">Liu et. al</a>] detailed survey of current estimates of the reproduction number. Parameters for the diseases' clinical characteristics are taken from the following <a href="https://www.who.int/docs/default-source/coronaviruse/who-china-joint-mission-on-covid-19-final-report.pdf">WHO Report</a>. 
</p>


<p class = "center">
<b> Model Details </b><br>
The clinical dynamics in this model are an elaboration on SEIR that simulates the disease's progression at a higher resolution, subdividing {@html math_inline("I,R")} into asymptomatic, symptomatic, hospitalized and patients in an ICU. Each of these variables follows its own trajectory to the final outcome, and the sum of these compartments add up to the values predicted by SEIR. Please refer to the source code for details. Note that we assume, for simplicity, that all fatalities come from ICU's, and that all fatal cases are admitted to hospitals.
</p>

<p class = "center">
<b> Acknowledgements </b><br>
<a href = "https://gabgoh.github.io/">Gabriel Goh</a> for the original version. 
<a href = "https://enkimute.github.io/">Steven De Keninck</a> for RK4 Integrator. <a href="https://twitter.com/ch402">Chris Olah</a>, <a href="https://twitter.com/shancarter">Shan Carter
</a> and <a href="https://twitter.com/ludwigschubert">Ludwig Schubert
</a> wonderful feedback. <a href="https://twitter.com/NikitaJer">Nikita Jerschov</a> for improving clarity of text. Charie Huang for context and discussion.
</p>

<!-- Input data -->
<div style="margin-bottom: 30px">

  <div class="center" style="padding: 10px; margin-top: 3px; width: 925px">
    <div class="legendtext">Export parameters:</div>
    <form>
      <textarea type="textarea" rows="1" cols="5000" style="white-space: nowrap;  overflow: auto; width:100%; text-align: left" id="fname" name="fname">{state}</textarea>
    </form>
  </div>
</div>
