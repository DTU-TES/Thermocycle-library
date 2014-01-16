within ThermoCycle.Components.HeatFlow.HeatTransfer.ConvectiveHeatTransfer;
model SmoothCorrelations
  "Liquid, two-phase and vapor correlations smoothened at phase transitions"
  extends
    ThermoCycle.Components.HeatFlow.HeatTransfer.ConvectiveHeatTransfer.BaseClasses.PartialConvectiveCorrelation;
  import SMOOTH = Modelica.Fluid.Dissipation.Utilities.Functions.General.Stepsmoother;

  parameter Modelica.SIunits.Length d_hyd = 0.010 "hydralic diameter";
  parameter Modelica.SIunits.Length L = 0.5 "cell length";
  parameter Modelica.SIunits.Length P = d_hyd*Modelica.Constants.pi
    "flow perimeter";
  parameter Modelica.SIunits.Area A_c = d_hyd^2/4*Modelica.Constants.pi
    "cross sectional area";

  Medium.SaturationProperties FluidSat[n];

  parameter Modelica.SIunits.CoefficientOfHeatTransfer U_init = 2000
    "initialization heat transfer coefficient";
  Modelica.SIunits.CoefficientOfHeatTransfer U(start = U_init)
    "heat transfer coefficient";

  parameter Modelica.SIunits.Time transitionStart = 5
    "start time of transition from nominal to correlation";
  parameter Modelica.SIunits.Time transitionTime = 5
    "transition time of transition from nominal to correlation";

  replaceable function HTC_liquid =
      ThermoCycle.Components.HeatFlow.HeatTransfer.ConvectiveHeatTransfer.BaseClasses.HeatTransferCorrelations.SinglePhaseCorrelations.Constant_htc
    constrainedby
    ThermoCycle.Components.HeatFlow.HeatTransfer.ConvectiveHeatTransfer.BaseClasses.HeatTransferCorrelations.SinglePhaseCorrelations.PartialSinglePhaseHeatTransfer(
    redeclare package Medium = Medium) annotation (choicesAllMatching=true);
  replaceable function HTC_vapor =
      ThermoCycle.Components.HeatFlow.HeatTransfer.ConvectiveHeatTransfer.BaseClasses.HeatTransferCorrelations.SinglePhaseCorrelations.Constant_htc
    constrainedby
    ThermoCycle.Components.HeatFlow.HeatTransfer.ConvectiveHeatTransfer.BaseClasses.HeatTransferCorrelations.SinglePhaseCorrelations.PartialSinglePhaseHeatTransfer(
    redeclare package Medium = Medium) annotation (choicesAllMatching=true);
  replaceable function HTC_twophase =
  ThermoCycle.Components.HeatFlow.HeatTransfer.ConvectiveHeatTransfer.BaseClasses.HeatTransferCorrelations.TwoPhaseCorrelations.Constant_htc
    constrainedby
    ThermoCycle.Components.HeatFlow.HeatTransfer.ConvectiveHeatTransfer.BaseClasses.HeatTransferCorrelations.TwoPhaseCorrelations.PartialTwoPhaseHeatTransfer(
    redeclare package Medium = Medium) annotation (choicesAllMatching=true);

protected
  Modelica.SIunits.CoefficientOfHeatTransfer Ucorr
    "vapor heat transfer coefficient";

equation
  FluidSat[1] = Medium.setSat_p(FluidState[1].p);
    // this could be passed by cell model, but computed here to not modify the cell models at this point..

  q_dot = {U*(thermalPortL[i].T - T_fluid[i]) for i in 1:n};

  if noEvent(time<transitionStart) then
     U = U_init;
     // I would have liked to use U = (Unom_l + Unom_tp + Unom_v)/3 as guess value in the declaration of U, ie.
     // U(start = (Unom_l + Unom_tp + Unom_v)/3)
     // However Unom is not passed as parameters, but becomes Real that cannot be used as guess value.
     // Could possibly eliminate U_init and use Unom here (and as guess value)...
  elseif noEvent(time>=transitionStart and time <= transitionStart+transitionTime) then
     U = SMOOTH(transitionStart,transitionStart+transitionTime,time)*U_init  + SMOOTH(transitionStart+transitionTime,transitionStart,time)*Ucorr;
  else
     U = Ucorr;
  end if;

  if noEvent(x<=0) then
    Ucorr = HTC_liquid(
      L=L,
      d_hyd=d_hyd,
      P=P,
      A_c=A_c,
      mdot=M_dot,
      state=FluidState[1]);
   elseif noEvent(x>0 and x<=0.05) then
     Ucorr = SMOOTH(
      0,
      0.05,
      x)*HTC_liquid(
      L=L,
      d_hyd=d_hyd,
      P=P,
      A_c=A_c,
      mdot=M_dot,
      state=Medium.setBubbleState(FluidSat[1])) + SMOOTH(
      0.05,
      0,
      x)*HTC_twophase(
      L=L,
      d_hyd=d_hyd,
      P=P,
      A_c=A_c,
      mdot=M_dot,
      q=sum(q_dot[i] for i in 1:n),
      x=x,
      sat=FluidSat[1]);
   elseif noEvent(x>0.05 and x<0.95) then
     Ucorr = HTC_twophase(
      L=L,
      d_hyd=d_hyd,
      P=P,
      A_c=A_c,
      mdot=M_dot,
      q=sum(q_dot[i] for i in 1:n),
      x=x,
      sat=FluidSat[1]);
   elseif noEvent(x>0.95 and x<1) then
     Ucorr = SMOOTH(
      1,
      0.95,
      x)*HTC_vapor(
      L=L,
      d_hyd=d_hyd,
      P=P,
      A_c=A_c,
      mdot=M_dot,
      state=Medium.setDewState(FluidSat[1])) + SMOOTH(
      0.95,
      1,
      x)*HTC_twophase(
      L=L,
      d_hyd=d_hyd,
      P=P,
      A_c=A_c,
      mdot=M_dot,
      q=sum(q_dot[i] for i in 1:n),
      x=x,
      sat=FluidSat[1]);
   else
     Ucorr = HTC_vapor(
      L=L,
      d_hyd=d_hyd,
      P=P,
      A_c=A_c,
      mdot=M_dot,
      state=FluidState[1]);
   end if;

  annotation(Documentation(info="<html>
<p>Simple heat transfer correlation with constant heat transfer coefficient. </p>
<p>Taken from: Modelica.Fluid.Pipes.BaseClasses.HeatTransfer.ConstantFlowHeatTransfer</p>
</html>"));
end SmoothCorrelations;
