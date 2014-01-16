within ThermoCycle.Examples.TestComponents;
model flow1D_R134a_with_corretions
  parameter Integer N = 10;
  ThermoCycle.Components.FluidFlow.Pipes.Flow1Dim flow1Dim(
    Unom_l=400,
    Unom_tp=1000,
    Unom_v=400,
    N=N,
    Discretization=ThermoCycle.Functions.Enumerations.Discretizations.upwind_AllowFlowReversal,
    redeclare package Medium = CoolProp2Modelica.Media.R134a_CP,
    A=0.628,
    V=0.00157,
    Mdotnom=0.03,
    pstart=1000000,
    Tstart_inlet=323.15,
    Tstart_outlet=373.15,
    redeclare model Flow1DimHeatTransferModel =
        ThermoCycle.Components.HeatFlow.HeatTransfer.ConvectiveHeatTransfer.SmoothCorrelations
        (
        redeclare function HTC_liquid =
            ThermoCycle.Components.HeatFlow.HeatTransfer.ConvectiveHeatTransfer.BaseClasses.HeatTransferCorrelations.SinglePhaseCorrelations.Gnielinski_MSL,

        redeclare function HTC_vapor =
            ThermoCycle.Components.HeatFlow.HeatTransfer.ConvectiveHeatTransfer.BaseClasses.HeatTransferCorrelations.SinglePhaseCorrelations.Gnielinski_MSL,

        L=20/10,
        redeclare function HTC_twophase =
            ThermoCycle.Components.HeatFlow.HeatTransfer.ConvectiveHeatTransfer.BaseClasses.HeatTransferCorrelations.TwoPhaseCorrelations.GungerWinterton1986_horizontalBoiling_MSL))
    annotation (Placement(transformation(extent={{-40,-4},{-2,34}})));

  ThermoCycle.Components.HeatFlow.Sources.Source_T source_T(N=N)
    annotation (Placement(transformation(extent={{-30,42},{4,64}})));
  Modelica.Blocks.Sources.Constant const(k=273.15 + 20)
    annotation (Placement(transformation(extent={{-62,64},{-54,72}})));
  Components.FluidFlow.Reservoirs.SourceMdot             sourceMdot1(
    UseT=false,
    redeclare package Medium = CoolProp2Modelica.Media.R134a_CP,
    h_0=1.9e5,
    Mdot_0=0.03,
    p=500000,
    T_0=356.26)
    annotation (Placement(transformation(extent={{-86,12},{-60,38}})));
  Components.FluidFlow.Reservoirs.SinkP             sinkP(
    h=254381,
    redeclare package Medium = CoolProp2Modelica.Media.R134a_CP,
    p0=500000)
    annotation (Placement(transformation(extent={{50,14},{70,34}})));
equation
  connect(source_T.thermalPort, flow1Dim.Wall_int) annotation (Line(
      points={{-13.17,48.49},{-13.17,43.95},{-21,43.95},{-21,22.9167}},
      color={255,0,0},
      smooth=Smooth.None));
  connect(const.y, source_T.Temperature) annotation (Line(
      points={{-53.6,68},{-13,68},{-13,57.4}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sourceMdot1.flangeB, flow1Dim.InFlow) annotation (Line(
      points={{-61.3,25},{-52,25},{-52,15},{-36.8333,15}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(flow1Dim.OutFlow, sinkP.flangeB) annotation (Line(
      points={{-5.16667,15.1583},{18,15.1583},{18,24},{51.6,24}},
      color={0,0,255},
      smooth=Smooth.None));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{
            100,100}}), graphics={Text(
          extent={{-62,56},{-26,50}},
          lineColor={0,0,0},
          textString="Thermal port")}),
    experiment(StopTime=50, __Dymola_Algorithm="Dassl"),
    __Dymola_experimentSetupOutput);
end flow1D_R134a_with_corretions;
