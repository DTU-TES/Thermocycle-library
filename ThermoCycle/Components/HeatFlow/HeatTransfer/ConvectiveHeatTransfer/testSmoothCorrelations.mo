within ThermoCycle.Components.HeatFlow.HeatTransfer.ConvectiveHeatTransfer;
model testSmoothCorrelations

 package Medium = ThermoCycle.Media.R134aCP;

 ThermoCycle.Interfaces.HeatTransfer.ThermalPortL Wall_int
    annotation (Placement(transformation(extent={{-28,40},{32,60}}),
        iconTransformation(extent={{-40,40},{40,60}})));

/************ Geometric characteristics **************/
  constant Real pi = Modelica.Constants.pi "pi-greco";
//  parameter Modelica.SIunits.Volume Vi "Volume of a single cell";
  parameter Modelica.SIunits.Area Ai=0.01*pi "Lateral surface of a single cell";
  parameter Modelica.SIunits.MassFlowRate Mdotnom=0.1 "Nominal fluid flow rate";
  parameter Modelica.SIunits.CoefficientOfHeatTransfer Unom_l=1500
    "if HTtype = LiqVap : Heat transfer coefficient, liquid zone ";
  parameter Modelica.SIunits.CoefficientOfHeatTransfer Unom_tp=3000
    "if HTtype = LiqVap : heat transfer coefficient, two-phase zone";
  parameter Modelica.SIunits.CoefficientOfHeatTransfer Unom_v=700
    "if HTtype = LiqVap : heat transfer coefficient, vapor zone";
 /************ FLUID INITIAL VALUES ***************/
  parameter Modelica.SIunits.Pressure pstart=5e5 "Fluid pressure start value"
                                     annotation (Dialog(tab="Initialization"));
  parameter Medium.SpecificEnthalpy hstart=2.1E5 "Start value of enthalpy"
    annotation (Dialog(tab="Initialization"));

replaceable model HeatTransfer =
      ThermoCycle.Components.HeatFlow.HeatTransfer.ConvectiveHeatTransfer.SmoothCorrelations
    constrainedby
    ThermoCycle.Components.HeatFlow.HeatTransfer.ConvectiveHeatTransfer.BaseClasses.PartialConvectiveCorrelation
    "Convective heat transfer"                                                         annotation (choicesAllMatching = true);
HeatTransfer heatTransfer(
redeclare final package Medium = Medium,
final n=1,
final Mdotnom = Mdotnom,
final Unom_l = Unom_l,
final Unom_tp = Unom_tp,
final Unom_v = Unom_v,
final M_dot = M_dot_su,
final x = x,
final FluidState={fluidState},
transitionStart = -5,
transitionTime = -5,
    redeclare function HTC_vapor =
        ThermoCycle.Components.HeatFlow.HeatTransfer.ConvectiveHeatTransfer.BaseClasses.HeatTransferCorrelations.SinglePhaseCorrelations.DittusBoelter_MSL,
    redeclare function HTC_twophase =
        ThermoCycle.Components.HeatFlow.HeatTransfer.ConvectiveHeatTransfer.BaseClasses.HeatTransferCorrelations.TwoPhaseCorrelations.GungerWinterton1986_horizontalBoiling_MSL,
    redeclare function HTC_liquid =
        ThermoCycle.Components.HeatFlow.HeatTransfer.ConvectiveHeatTransfer.BaseClasses.HeatTransferCorrelations.SinglePhaseCorrelations.Gnielinski_MSL,
    U_init=2000)       annotation (Placement(transformation(extent={{-12,-14},{8,6}})));

/***************  VARIABLES ******************/
  Medium.ThermodynamicState  fluidState;

 // Medium.ThermodynamicState State1;
  Medium.SaturationProperties sat;
  Medium.AbsolutePressure p(start=pstart);
  Medium.SpecificEnthalpy h(start=hstart)
    "Fluid specific enthalpy at the cells";
  Medium.Temperature T "Fluid temperature";
  Medium.Density rho "Fluid cell density";
  Modelica.SIunits.HeatFlux qdot "heat flux at each cell";
  Real x "Vapor quality";

  Modelica.SIunits.MassFlowRate M_dot_su;
  Modelica.SIunits.HeatFlowRate Q_tot;
  Modelica.SIunits.SpecificEnthalpy h_v;
  Modelica.SIunits.SpecificEnthalpy h_l;

equation
  M_dot_su = Mdotnom;
  der(p)=0;
  der(h)=2.5e5;

  sat = Medium.setSat_p(p);
  h_v = Medium.dewEnthalpy(sat);
  h_l = Medium.bubbleEnthalpy(sat);
  //T_sat = Medium.temperature(sat);
  /* Fluid Properties */
  fluidState = Medium.setState_ph(p,h);
  T = Medium.temperature(fluidState);
  rho = Medium.density(fluidState);
  x = (h - h_l)/(h_v - h_l);
  qdot = heatTransfer.q_dot[1];
  Q_tot = Ai*qdot;
  connect(heatTransfer.thermalPortL[1], Wall_int) annotation (Line(
      points={{-2.2,2.6},{-2.2,28.3},{2,28.3},{2,50}},
      color={255,0,0},smooth=Smooth.None));

  annotation (Diagram(graphics));
end testSmoothCorrelations;
