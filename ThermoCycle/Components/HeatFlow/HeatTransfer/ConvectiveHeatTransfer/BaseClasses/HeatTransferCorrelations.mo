within ThermoCycle.Components.HeatFlow.HeatTransfer.ConvectiveHeatTransfer.BaseClasses;
package HeatTransferCorrelations
  extends Modelica.Icons.Package;

  package SinglePhaseCorrelations
      extends Modelica.Icons.Package;

    partial function PartialSinglePhaseHeatTransfer
      "Partial function with main inputs and output"

      replaceable package Medium =
          Modelica.Media.Interfaces.PartialMedium "Medium model";

      input Modelica.SIunits.Length L "flow length";
      input Modelica.SIunits.Length d_hyd "hydraulic diameter";
      input Modelica.SIunits.Length P "perimeter";
      input Modelica.SIunits.Area A_c "cross area";

      input Modelica.SIunits.MassFlowRate mdot;
      input Medium.ThermodynamicState state;

      output Modelica.SIunits.CoefficientOfHeatTransfer htc;

    end PartialSinglePhaseHeatTransfer;

    function DittusBoelter_MSL
      "Dittus Boelter correlation from Modelica Standard Library - Only internal pipe flow"

      extends PartialSinglePhaseHeatTransfer;

    protected
    Modelica.Fluid.Dissipation.HeatTransfer.StraightPipe.kc_overall_IN_con IN_con(
        target=Modelica.Fluid.Dissipation.Utilities.Types.HeatTransferBoundary.UWTuDFF,
        d_hyd=d_hyd,
        L=L,
        roughness=Modelica.Fluid.Dissipation.Utilities.Types.Roughness.Neglected)
        annotation (Placement(transformation(extent={{52,60},{72,80}})));

    Modelica.Fluid.Dissipation.HeatTransfer.StraightPipe.kc_overall_IN_var
        IN_var(
        cp=Medium.specificHeatCapacityCp(state),
        eta=Medium.dynamicViscosity(state),
        rho=Medium.density(state),
        lambda=Medium.thermalConductivity(state),
        m_flow=mdot)
               annotation (Placement(transformation(extent={{20,60},{40,80}})));

    algorithm
    htc := Modelica.Fluid.Dissipation.HeatTransfer.StraightPipe.kc_overall_KC(IN_con,IN_var);

    end DittusBoelter_MSL;

    function Gnielinski_MSL
      "Gnielinski correlation from Modelica Standard Library - Only internal pipe flow"

      extends PartialSinglePhaseHeatTransfer;

    protected
    Modelica.Fluid.Dissipation.HeatTransfer.StraightPipe.kc_overall_IN_con IN_con(
        target=Modelica.Fluid.Dissipation.Utilities.Types.HeatTransferBoundary.UWTuDFF,
        d_hyd=d_hyd,
        L=L,
        roughness=Modelica.Fluid.Dissipation.Utilities.Types.Roughness.Considered)
        annotation (Placement(transformation(extent={{52,60},{72,80}})));

    Modelica.Fluid.Dissipation.HeatTransfer.StraightPipe.kc_overall_IN_var
        IN_var(
        cp=Medium.specificHeatCapacityCp(state),
        eta=Medium.dynamicViscosity(state),
        rho=Medium.density(state),
        lambda=Medium.thermalConductivity(state),
        m_flow=mdot)
               annotation (Placement(transformation(extent={{20,60},{40,80}})));

    algorithm
    htc := Modelica.Fluid.Dissipation.HeatTransfer.StraightPipe.kc_overall_KC(IN_con,IN_var);

    end Gnielinski_MSL;

    function Constant_htc "Constant heat transfer coefficient"

      extends PartialSinglePhaseHeatTransfer;

      parameter Modelica.SIunits.CoefficientOfHeatTransfer htc_const=800;

    algorithm
      htc := htc_const;

    end Constant_htc;
  end SinglePhaseCorrelations;

  package TwoPhaseCorrelations
      extends Modelica.Icons.Package;

    partial function PartialTwoPhaseHeatTransfer
      replaceable package Medium =
          Modelica.Media.Interfaces.PartialTwoPhaseMedium "Medium model";

      input Modelica.SIunits.Length L "flow length";
      input Modelica.SIunits.Length d_hyd "hydraulic diameter";
      input Modelica.SIunits.Length P "perimeter";
      input Modelica.SIunits.Area A_c "cross area";

      input Modelica.SIunits.MassFlowRate mdot;
      input Modelica.SIunits.HeatFlux q;
      input Modelica.SIunits.QualityFactor x;
      input Medium.SaturationProperties sat;

      output Modelica.SIunits.CoefficientOfHeatTransfer htc;

    end PartialTwoPhaseHeatTransfer;

    function GungerWinterton1986_horizontalBoiling_MSL
      "Gunger and Winterton 1986 (horizontal flow boiling) from Modelica Standard Library"

      extends PartialTwoPhaseHeatTransfer;

    protected
    Modelica.Fluid.Dissipation.HeatTransfer.StraightPipe.kc_twoPhaseOverall_KC_IN_con
                                                                                        IN_con(
        A_cross=A_c,
        perimeter=P,
        p_crit=Medium.fluidConstants[1].criticalPressure,
        MM=Medium.fluidConstants[1].molarMass)                                                                                                     annotation (Placement(transformation(extent={{52,60},{72,80}})));
    // not correct for mixtures..

    Medium.ThermodynamicState state_l = Medium.setBubbleState(sat);
    Medium.ThermodynamicState state_v = Medium.setDewState(sat);

    Modelica.Fluid.Dissipation.HeatTransfer.StraightPipe.kc_twoPhaseOverall_KC_IN_var
                                                                                      IN_var(
        cp_l=Medium.specificHeatCapacityCp(state_l),
        lambda_l=Medium.thermalConductivity(state_l),
        rho_l=Medium.bubbleDensity(sat),
        rho_g=Medium.dewDensity(sat),
        eta_l=Medium.dynamicViscosity(state_l),
        eta_g=Medium.dynamicViscosity(state_v),
        dh_lg=Medium.dewEnthalpy(sat)-Medium.bubbleEnthalpy(sat),
        m_flow=mdot,
        qdot_A=q,
        x_flow=x,
        pressure=sat.psat) annotation (Placement(transformation(extent={{20,60},{40,80}})));

    algorithm
        htc :=
        Modelica.Fluid.Dissipation.Utilities.Functions.HeatTransfer.TwoPhase.kc_twoPhase_boilingHorizontal_KC(
        IN_con, IN_var);

    end GungerWinterton1986_horizontalBoiling_MSL;

    function GungerWinterton1986_verticalBoiling_MSL
      "Gunger and Winterton 1986 (vertical flow boiling) from Modelica Standard Library"

      extends PartialTwoPhaseHeatTransfer;

    protected
    Modelica.Fluid.Dissipation.HeatTransfer.StraightPipe.kc_twoPhaseOverall_KC_IN_con
                                                                                        IN_con(
        A_cross=A_c,
        perimeter=P,
        p_crit=Medium.fluidConstants[1].criticalPressure,
        MM=Medium.fluidConstants[1].molarMass)                                                                                                     annotation (Placement(transformation(extent={{52,60},{72,80}})));
    // not correct for mixtures..

    Medium.ThermodynamicState state_l = Medium.setBubbleState(sat);
    Medium.ThermodynamicState state_v = Medium.setDewState(sat);

    Modelica.Fluid.Dissipation.HeatTransfer.StraightPipe.kc_twoPhaseOverall_KC_IN_var
                                                                                      IN_var(
        cp_l=Medium.specificHeatCapacityCp(state_l),
        lambda_l=Medium.thermalConductivity(state_l),
        rho_l=Medium.bubbleDensity(sat),
        rho_g=Medium.dewDensity(sat),
        eta_l=Medium.dynamicViscosity(state_l),
        eta_g=Medium.dynamicViscosity(state_v),
        dh_lg=Medium.dewEnthalpy(sat)-Medium.bubbleEnthalpy(sat),
        m_flow=mdot,
        qdot_A=q,
        x_flow=x,
        pressure=sat.psat) annotation (Placement(transformation(extent={{20,60},{40,80}})));

    algorithm
        htc :=
        Modelica.Fluid.Dissipation.Utilities.Functions.HeatTransfer.TwoPhase.kc_twoPhase_boilingVertical_KC(
        IN_con, IN_var);

    end GungerWinterton1986_verticalBoiling_MSL;

    function Shah1979_horizontalCondensation_MSL
      "Shah 1979 (horizontal flow condensation) from Modelica Standard Library"

      extends PartialTwoPhaseHeatTransfer;

    protected
    Modelica.Fluid.Dissipation.HeatTransfer.StraightPipe.kc_twoPhaseOverall_KC_IN_con
                                                                                        IN_con(
        A_cross=A_c,
        perimeter=P,
        p_crit=Medium.fluidConstants[1].criticalPressure,
        MM=Medium.fluidConstants[1].molarMass)                                                                                                     annotation (Placement(transformation(extent={{52,60},{72,80}})));
    // not correct for mixtures..

    Medium.ThermodynamicState state_l = Medium.setBubbleState(sat);
    Medium.ThermodynamicState state_v = Medium.setDewState(sat);

    Modelica.Fluid.Dissipation.HeatTransfer.StraightPipe.kc_twoPhaseOverall_KC_IN_var
                                                                                      IN_var(
        cp_l=Medium.specificHeatCapacityCp(state_l),
        lambda_l=Medium.thermalConductivity(state_l),
        rho_l=Medium.bubbleDensity(sat),
        rho_g=Medium.dewDensity(sat),
        eta_l=Medium.dynamicViscosity(state_l),
        eta_g=Medium.dynamicViscosity(state_v),
        dh_lg=Medium.dewEnthalpy(sat)-Medium.bubbleEnthalpy(sat),
        m_flow=mdot,
        qdot_A=q,
        x_flow=x,
        pressure=sat.psat) annotation (Placement(transformation(extent={{20,60},{40,80}})));

    algorithm
        htc :=
        Modelica.Fluid.Dissipation.Utilities.Functions.HeatTransfer.TwoPhase.kc_twoPhase_condensationHorizontal_KC(
        IN_con, IN_var);

    end Shah1979_horizontalCondensation_MSL;

    function Constant_htc "Constant heat transfer correlation"

      extends PartialTwoPhaseHeatTransfer;

      parameter Modelica.SIunits.CoefficientOfHeatTransfer htc_const=3000;

    algorithm
      htc := htc_const;

    end Constant_htc;
  end TwoPhaseCorrelations;


end HeatTransferCorrelations;
