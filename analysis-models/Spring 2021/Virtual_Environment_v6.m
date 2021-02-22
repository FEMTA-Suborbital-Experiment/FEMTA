clear
clc
close all
load('altitude.mat');
load('time.mat');

%CONSTANT AND INITIAL VALUES
global R kB N_a;
R=8.3145;                                     %Universal Gas Constant [J/mol-K]
kB=1.380649e-23;                              %Boltzmann constant [J/K]
N_a=6.02214e23;                               %Avogadros number [# of particles/mol]

%Air Properties
MW_Air=29;                                    %Molecular Weight of Air [g/mol]
R_air=R/MW_Air;                               %Specific Gas Constant for Air [J/g-K]
Cp_Air=1.005;                                 %Specific Heat of Air [kJ/kg-K]
gammaAir=Cp_Air/(Cp_Air-R_air);               %ratio of specific heats for air

%Water Properties
rho_water=997;                                %density of liquid water [kg/m^3]
m_H2O=2.988e-26;                              %mass of one water molecule [kg]
MW_Water=18.0135;                             %Molecular weight of water [g/mol]
R_wv=R/MW_Water;                              %Specific Gas Constant for Water Vapor [J/g-K]
Cp_water_liquid=4.2174;                       %specific heat of liquid water [kJ/kg-K]
Cp_water_vapor=1.864;                         %specific heat of water vapor [kJ/kg-K]
Ce_water=1;                                   %Evaporation Coefficient of water
Cc_water=0.75;                                %Condensation Coefficient of water
gammaWV=Cp_water_vapor/(Cp_water_vapor-R_wv); %ratio of specific heats for water vapor

%HFE Properties
MW_HFE=250;                                   %Molecular Weight of HFE [g/mol]
R_HFE=R/MW_HFE;                               %Specific Gas Constant for HFE [J/g-K]
h_evap_HFE=125.6;                             %Heat of Vaporization of HFE [kJ/kg]
Cp_HFEliquid=1.172303;                        %Specific Heat of HFE liquid [kJ/kg-K]
m_HFE=MW_HFE/(N_a*1000);                      %Mass of one HFE molecule [kg]
Ce_HFE=1;                                     %Evaporation Coefficient of HFE
Cc_HFE=0.1;                                   %Condensation Coefficient of HFE

%Propellant Tank
P0_tank=101325;                                         %Initial pressure in prop tank [Pa] (1 atm)
T0_tank=300;                                            %Initial temperature in prop tank [K]
V_tank=142.567e-6;                                      %Total volume of both prop tanks [m^3] (142.567mL) (CONSTANT)
volAir0=0.328e-6;                                       %Initial volume of air in both prop tanks [m^3] (0.328)
nAir_tank=(P0_tank*volAir0)/(R*T0_tank);                %Number of moles of air in both prop tanks (CONSTANT)
volHFE_liquid0_tank=0.983e-6;                           %initial volume of HFE in both prop tanks [m^3] (0.983mL)
volWater0_tank=V_tank-volAir0-volHFE_liquid0_tank;      %initial volume of water in both prop tanks [m^3]
A_HFE=3.167e-5;                                         %Area from which HFE condenses and evaporates [m^2]

%Collection Chamber
P0_CC=101325;                                           %initial pressure in collection chamber [Pa]
T0_CC=300;                                              %initial temperature in collection chamber [K]
V_CC=334.519e-6;                                        %volume of CC [m^3] (542.248mL)
ventSolenoidDiam=2.18e-3;                               %Daimeter of vent solenoid [m]
CCBeta=ventSolenoidDiam/0.09398;                        %Ratio of vent solenoid orifice to CC cross section
nAir_CC=(P0_CC*V_CC)/(R*T0_CC);                         %Initial number of moles of air in CC

%Piping network
D_pipe=(1/8)/39.37;       %Pipe diameter [m] (CONSTANT)
A=pi*((D_pipe/2)^2);      %Area of pipe cross section [m^2] (CONSTANT)

%Orifice
D_O=0.5e-3;               %Orifice diameter [m]
A_O=pi*((D_O/2).^2);      %Area of orifice [m^2]
Beta=D_O/D_pipe;          %Ratio of orifice to pipe diameter
CD_orifice=0.6;           %discharge coefficient of orifice

%VARIABLES
volWater_tank=volWater0_tank;                      %Initial volume of water in one Prop Tank [m^3]
volWater_CC=0;                                     %Initial volume of water in Collection Chamber [m^3]
volHFE_liquid=volHFE_liquid0_tank;                 %Volume of HFE Liquid in Prop Tank [m^3]
tankPress=P0_tank;                                 %Pressure in Prop Tank [Pa]
CCPress=P0_CC;                                     %Pressure in Collection Chamber [Pa]
tankTempGas=T0_tank;                               %Temperature of gas in Prop Tank [K]
tankTempLiquid_HFE=T0_tank;                        %Temperature of HFE Liquid in Prop Tank [K]
CCTempGas=T0_CC;                                   %Temperature of gas in CC [K]
CCTempLiquid=T0_CC;                                %Temperature of Water Liquid in CC [K]
m_HFE_vapor=0;                                     %Initial mass of HFE vapor in prop tank
m_HFE_liquid=volHFE_liquid0_tank*nvcRho(T0_tank);  %Initial mass of HFE liquid in prop tank
m_water_vapor=0;                                   %Initial mass of water vapor in CC
m_water_liquid=0;                                  %Initial mass of water liquid in CC
n_Gas=nAir_tank;                                   %Initial moles of gas in tank
n_HFE_vapor=0;                                     %Initial moles of HFE vapor in tank
volGas=volAir0;                                    %Initial volume of gas in tank [m^3]
volWater_shut=0;                                   %Initial volume of water in prop tank when valve shuts [m^3]
n_air_lost=0;
A_HFE_cond=A_HFE;
A_HFE_evap=A_HFE;
time=0;
count=1;                           
dt=2e-4;    %timestep [s] (do not run at more than 2e-4) (loop time: 255s @ 2e-4)

%ARRAY INITIALIZATION
tankVolWater_array=[];
tankPress_array=[];
tankTempGas_array=[];
tankTempLiquid_array=[];
CCVolWater_array=[];
CCPress_array=[];
CCPress_preExp_array=[];
CCTempGas_array=[];
CCTempLiquid_array=[];
PvapHFE_array=[];
PvapWater_array=[];
QHFE_array=[];
Qwater_array=[];
time_array=[];
m_HFE_transfer_array=[];
m_HFE_vapor_array=[];
m_HFE_liquid_array=[];
tankVolGas_array=[];
n_Gas_array=[];
m_water_vapor_array=[];
m_water_liquid_array=[];
m_water_transfer_array=[];
m_water_unaltered_array=[];
h_evap_water_array=[];
m_HFE_unaltered_array=[];
volGas_array=[];
nGas_array=[];
m_HFE_total_array=[];
m_water_total_array=[];
m_water_liquid_tank_array=[];
nAir_CC_array=[];
m_water_lost_array=[];
n_Air_lost_array=[];
m_lost_array=[];
nWaterVapor_CC_array=[];
altitude_array=[];
flo_water_array=[];
exp_time_array=[];

%EXECUTABLE LOOPS
tic
while time<max(t)
        
        alt=interp1(t,h,time);
        [ambientT,ambientP,ambientRho]=StandardAtm(alt);
        
        if volWater_tank-volWater_shut<0
            fprintf("Tank Volume Empty\n\n")
            break;
        end
        
        %Condition for beginning experiment
        if alt<80000
            ventSol=1;
            flowSol=0;
        else
            flowSol=2;
            ventSol=1;
            if time>200
                flowSol=1;
            end
        end
        
        %Conditions for loop termination
        if flowSol==1 && volWater_shut==0
            volWater_shut=0.5*volWater_tank;
        elseif flowSol==2
            volWater_shut=0;
        end

        %Volumetric Flow Rate of Liquid Water Propellant through one orifice[m^3/s]
        flo_water=flowSol*CD_orifice*A_O*sqrt(2*(max(tankPress,CCPress)-min(tankPress,CCPress))/(rho_water*(1-(Beta^4))));
        if tankPress<CCPress
            flo_water=flo_water*-1;
        end

        
        if ~isreal(flo_water) || isnan(flo_water)
            flo_water
            fprintf("\nError: flow_water imaginary or NaN\n\n")
            break;
        end

        %Update of Liquid Water Propellant Volumes in Tank & CC [m^3]
        volWater_tank=volWater_tank-(flo_water*dt);
        volWater_CC=volWater_CC+(flo_water*dt);
        m_water_liquid_tank=volWater_tank*rho_water;

        %PROPELLANT TANK
            %Mass of Gas in the tank [kg]
            mAir_tank=(nAir_tank*MW_Air)/1000;
            xAir2_tank=nAir_tank/(nAir_tank+n_HFE_vapor);
            xHFEvapor2_tank=n_HFE_vapor/(nAir_tank+n_HFE_vapor);
            yAir2_tank=mAir_tank/(mAir_tank+m_HFE_vapor);
            yHFEvapor2_tank=m_HFE_vapor/(mAir_tank+m_HFE_vapor);
            MWGas2_tank=(MW_Air*xAir2_tank)+(MW_HFE*xHFEvapor2_tank);
            mGas2_tank=(n_Gas*MWGas2_tank)/1000;
            
            %Specific Heat of Gas in Tank before HFE transfer (state 2) [kJ/kg-K]
            CpGas2_tank=(Cp_Air*yAir2_tank)+(Cp_HFEliquid*yHFEvapor2_tank);
            
            %Vapor Pressure of HFE [Pa]
            Pvap_HFE=nvcVP(tankTempLiquid_HFE);

            %Density of HFE liquid [kg/m^3]
            rho_HFE=nvcRho(tankTempLiquid_HFE);

            %Pressure of Gas [Pa]
            tankPress=(n_Gas*R*tankTempGas)/volGas;

            %Amount of HFE either condensing or evaporating [kg]
            if m_HFE_liquid==0
                A_HFE_evap=0;
            elseif m_HFE_vapor==0
                A_HFE_cond=0;
            end
            m_HFE_transfer=HerKnu(Pvap_HFE,tankTempLiquid_HFE,tankTempGas,tankPress,m_HFE,A_HFE_evap,A_HFE_cond,Ce_HFE,Cc_HFE)*dt;
            if tankPress>Pvap_HFE && m_HFE_transfer>0
                m_HFE_transfer=0;
            end
            
            %Update mass of HFE liquid and vapor [kg]
            m_HFE_liquid=m_HFE_liquid-m_HFE_transfer;
            m_HFE_vapor=m_HFE_transfer+m_HFE_vapor;
            
            %Update moles/volume of HFE vapor/liquid 
            vol_HFE_liquid=m_HFE_liquid/rho_HFE;
            n_HFE_vapor=(m_HFE_vapor*1000)/MW_HFE;
            
            %Temperature of Gas [K]
            tankTempGas=((m_HFE_transfer*Cp_HFEliquid*tankTempLiquid_HFE)+(mGas2_tank*CpGas2_tank*tankTempGas))/((Cp_HFEliquid*m_HFE_transfer)+(mGas2_tank*CpGas2_tank));

            %Update total amount of Gas (Air + HFE)
            n_Gas=n_HFE_vapor+nAir_tank;
            volGas=V_tank-volWater_tank-vol_HFE_liquid;

            %Temperatrue Update [K]
            Q_HFE=m_HFE_transfer*h_evap_HFE;
            tankTempLiquid_HFE=(-Q_HFE/(m_HFE_liquid*Cp_HFEliquid))+tankTempLiquid_HFE;

            %Check total mass of HFE in Prop Tank [kg]
            m_HFE_total=m_HFE_liquid+m_HFE_vapor;

        %COLLECTION CHAMBER
            %Surface Area of Collected Water [m^2]
            r=((3*volWater_CC)/(4*pi))^(1/3);
            A_water=4*pi*(r^2);

            %Vapor Pressure of Water [Pa]
            Pvap_water=waterVP(CCTempLiquid);

            %Evaporation Heat of Water [kJ/kg]
            h_evap_water=waterHV(CCTempLiquid)/1000;

            %Mass of water either evaporating or condensing at current timestep [kg]
            m_water_transfer=HerKnu(Pvap_water,CCTempLiquid,CCTempGas,CCPress,m_H2O,A_water,A_water,Ce_water,Cc_water)*dt;
            if CCPress>Pvap_water && m_water_transfer>0
                m_water_transfer=0;
            end
            
            %Mass of gas lost through vent solenoid [kg]
            nWaterVapor_CC=(m_water_vapor*1000)/MW_Water;
            rhoGas_CC=(((nAir_CC*MW_Air)/1000)+((nWaterVapor_CC*MW_Water)/1000))/(V_CC-(m_water_liquid/rho_water));
            xAir_CC=nAir_CC/(nAir_CC+nWaterVapor_CC);
            xWaterVapor_CC=nWaterVapor_CC/(nWaterVapor_CC+nAir_CC);
            gammaGas_CC=1+(1/((xWaterVapor_CC/(gammaWV-1))+(xAir_CC/(gammaAir-1))));
            m_lost=mDotThruOrifice(CCPress,ambientP,rhoGas_CC,gammaGas_CC,0.1,ventSolenoidDiam*ventSol)*dt;
            m_water_lost=m_lost*xWaterVapor_CC;
            n_Air_lost=((m_lost*1000)*xAir_CC)/MW_Air;
            
            %Moles of Air in CC 
            nAir_CC=nAir_CC-n_Air_lost;

            %Total mass of water vapor and liquid at current time [kg]
            m_water_vapor=m_water_transfer-m_water_lost+m_water_vapor;
            m_water_liquid=(volWater_CC*rho_water)-m_water_vapor;

            %Pressure Update [Pa]
            CCPress=((nAir_CC+nWaterVapor_CC)*R*CCTempGas)/(V_CC-(m_water_liquid/rho_water));
            
            %Liquid Temperature Update [K]
            if m_water_liquid~=0
                Q_water=m_water_transfer*h_evap_water;
                mwn=flo_water*rho_water*dt;
                T1=((m_water_liquid-mwn)/m_water_liquid)*CCTempLiquid;
                T2=(mwn/m_water_liquid)*300;
                T3=-Q_water/(m_water_liquid*Cp_water_liquid);
                CCTempLiquid=T1+T2+T3;
            end

        %Update total mass of water in the system [kg]
        m_water_total=m_water_vapor+m_water_liquid+(volWater_tank*rho_water);

        %ARRAY UPDATE
            %Propellant Tank
            tankVolWater_array(count)=volWater_tank;
            tankPress_array(count)=tankPress;
            tankTempGas_array(count)=tankTempGas;
            tankTempLiquid_array(count)=tankTempLiquid_HFE;
            tankVolGas_array(count)=volGas;
            m_HFE_vapor_array(count)=m_HFE_vapor;
            m_HFE_liquid_array(count)=m_HFE_liquid;
            n_Gas_array(count)=n_Gas;
            volGas_array(count)=volGas;
            nGas_array(count)=n_Gas;
            m_HFE_total_array(count)=m_HFE_total;
            m_water_liquid_tank_array(count)=m_water_liquid_tank;
            m_HFE_transfer_array(count)=m_HFE_transfer;

            %Collection Chamber
            CCVolWater_array(count)=volWater_CC;
            CCPress_array(count)=CCPress;
            CCTempGas_array(count)=CCTempGas;
            CCTempLiquid_array(count)=CCTempLiquid;
            m_water_vapor_array(count)=m_water_vapor;
            m_water_liquid_array(count)=m_water_liquid;
            m_water_transfer_array(count)=m_water_transfer;
            h_evap_water_array(count)=h_evap_water;
            m_water_total_array(count)=m_water_total;
            nAir_CC_array(count)=nAir_CC;
            m_water_lost_array(count)=m_water_lost;
            n_Air_lost_array(count)=n_Air_lost;
            m_lost_array(count)=m_lost;
            nWaterVapor_CC_array(count)=nWaterVapor_CC;

            %Miscellaneous
            PvapHFE_array(count)=Pvap_HFE;
            QHFE_array(count)=Q_HFE;
            PvapWater_array(count)=Pvap_water;
            flo_water_array(count)=flo_water;
            time_array(count)=time;
            exp_time_array(count)=exp_time;
            altitude_array(count)=alt;

            time=time+dt;
            count=count+1;
end
toc

%Array Transpose
time_array=time_array';
exp_time_array=exp_time_array';
altitude_array=altitude_array';
m_water_total_array=m_water_total_array';
m_HFE_total_array=m_HFE_total_array';
volGas_array=volGas_array';
nGas_array=nGas_array';
tankVolWater_array=tankVolWater_array';
tankPress_array=tankPress_array';
tankTempGas_array=tankTempGas_array';
tankTempLiquid_array=tankTempLiquid_array';
tankVolGas_array=tankVolGas_array';
m_HFE_vapor_array=m_HFE_vapor_array';
m_HFE_liquid_array=m_HFE_liquid_array';
n_Gas_array=n_Gas_array';
m_HFE_unaltered_array=m_HFE_unaltered_array';
CCVolWater_array=CCVolWater_array';
CCPress_array=CCPress_array';
CCTempGas_array=CCTempGas_array';
CCTempLiquid_array=CCTempLiquid_array';
m_water_vapor_array=m_water_vapor_array';
m_water_liquid_array=m_water_liquid_array';
m_water_transfer_array=m_water_transfer_array';
m_water_unaltered_array=m_water_unaltered_array';
h_evap_water_array=h_evap_water_array';
PvapHFE_array=PvapHFE_array';
QHFE_array=QHFE_array';
PvapWater_array=PvapWater_array';
Qwater_array=Qwater_array';
m_HFE_transfer_array=m_HFE_transfer_array';
CCPress_preExp_array=CCPress_preExp_array';
m_water_liquid_tank_array=m_water_liquid_tank_array';
nAir_CC_array=nAir_CC_array';
m_water_lost_array=m_water_lost_array';
m_lost_array=m_lost_array';
nWaterVapor_CC_array=nWaterVapor_CC_array';
flo_water_array=flo_water_array';
n_Air_lost_array=n_Air_lost_array';
flo_water_exp_array=nonzeros(flo_water_array);

%PLOTS
%Volumes of both Prop Tank and CC
figure(1)
plot(time_array,tankVolWater_array*10^6,'Linewidth',3)
hold on
plot(time_array,CCVolWater_array*10^6,'Linewidth',3)
xlabel("Time [s]",'Fontsize',17)
ylabel("Volume of Water [mL]",'Fontsize',17)
legend("Propellant Tank","Collection Chamber",'Fontsize',15)
title("Volumes of Collection Chamber and Propellant Tank",'Fontsize',22)

%Temperature in Prop Tank
figure(2)
plot(time_array,tankTempGas_array,'Linewidth',3)
hold on
plot(time_array,tankTempLiquid_array,'Linewidth',3)
xlabel("Time [s]",'Fontsize',17)
ylabel("Temperature [K]",'Fontsize',17)
legend("Gas (Air + HFE Vapor)","Liquid HFE",'Fontsize',15)
title("Temperature in Propellant Tank",'Fontsize',22)

%Temperature in CC
figure(3)
plot(time_array,CCTempGas_array,'Linewidth',3)
hold on
plot(time_array,CCTempLiquid_array,'Linewidth',3)
legend("Water Vapor","Liquid Water",'Fontsize',15)
title("Temperature in Collection Chamber",'Fontsize',22)
xlabel("Time [s]",'Fontsize',17)
ylabel("Temperature [K]",'Fontsize',17)

%Pressure in Prop Tank
figure(4)
plot(time_array,tankPress_array/1000,'Linewidth',3)
hold on
plot(time_array,PvapHFE_array/1000,'Linewidth',3)
xlabel("Time [s]",'Fontsize',17)
ylabel("Pressure [kPa]",'Fontsize',17)
legend("Propellant Tank Pressure","Vapor Pressure of HFE",'Fontsize',15)
title("Pressures in Propellant Tank",'Fontsize',22)

%Pressures in CC
figure(5)
plot(time_array,CCPress_array/1000,'Linewidth',3)
hold on
plot(time_array,PvapWater_array/1000,'Linewidth',3)
xlabel("Time [s]",'Fontsize',17)
ylabel("Pressure [kPa]",'Fontsize',17)
legend("Collection Chamber Pressure","Vapor Pressure of Water",'Fontsize',15)
title("Pressures in Collection Chamber",'Fontsize',22)

%Mass of HFE in All States
figure(6)
plot(time_array,m_HFE_vapor_array*1000,'Linewidth',2)
hold on
plot(time_array,m_HFE_liquid_array*1000,'Linewidth',2)
plot(time_array,m_HFE_total_array*1000,'Linewidth',2)
xlabel("Time [s]",'Fontsize',17)
ylabel("Mass of HFE [g]",'Fontsize',17)
title("Amount of HFE in Each State",'Fontsize',22)
legend("HFE Vapor","HFE Liquid","HFE Total",'Fontsize',15)

%Mass of Water in All States
figure(7)
plot(time_array,m_water_vapor_array*1000,'Linewidth',2)
hold on
plot(time_array,m_water_liquid_array*1000,'Linewidth',2)
plot(time_array,m_water_total_array*1000,'Linewidth',2)
plot(time_array,tankVolWater_array*rho_water*1000,'Linewidth',2)
xlabel("Time [s]",'Fontsize',17)
ylabel("Mass of Water [g]",'Fontsize',17)
title("Amount of Water in Each State",'Fontsize',22)
legend("Water Vapor","Water Liquid","Water Total","Water in One Prop Tank",'Fontsize',15)

%Altitude vs Time
figure(8)
plot(time_array,altitude_array,'Linewidth',2)
xlabel("Time [s]",'Fontsize',17)
ylabel("Altitude [m]",'Fontsize',17)
title("Flight Altitude",'Fontsize',22)

%Vapor Pressure of NV 7100 (T in K, vp in Pa)
function vp = nvcVP(T)
    vp = exp(22.415 - 3641.9 * (1/T));
end

%Density of NV 7100 Liquid (T in K)
function dD = nvcRho(T)
    dD = 1.5383 - 0.002269*(T-273.15);
    dD = dD / (1000 * 0.000001); %kg/m^3
end

%Vapor Pressure of Water (T in K, vp in Pa)
function vp = waterVP(T)
    vp=10.^(8.07131-(1730.63./(233.426+(T-273))));        %pressure in mmHg
    vp=vp*133;                                            %conversion to Pa from mmHg
end

%Mass transfer from gas to liquid [kg/s] (HERTZ-KNUDSEN EQUATION)
%Negative denotes vapor to liquid (condensation), Positive denotes liquid to vapor
%(evaporation)
function m_transfer=HerKnu(Ps,T_liquid,T_vapor,Pg,m,A_evap,A_cond,C_evap,C_cond)
    global kB;
    m_transfer=(sqrt(m/(2*pi*kB))*((A_evap*C_evap*(Ps/sqrt(T_liquid)))-(A_cond*C_cond*(Pg/sqrt(T_vapor)))));
end

%Heat of Vaporization of Water [J/kg]
function hv=waterHV(T)
    Hvs=[2500.9 2496.2 2491.4 2477.2 2467.7 2458.3 2453.5 2441.7 2429.8 2420.3 2406 2396.4 2381.9 2372.3 2357.7 2333 2308 2282.5 2266.9 2256.4 2229.6 2202.1 2144.3 2082 2014.2 1939.7 1857.4 1765.4 1661.6 1543 1404.6 1238.4 1027.3 719.8];
    Ts=[0.00 2 4 10 14 18 20 25 30 34 40 44 50 54 60 70 80 90 96 100 110 120 140 160 180 200 220 240 260 280 300 320 340 360]+273;
    hv=interp1(Ts, Hvs, T, 'linear').*1000;
    if isnan(hv)
        hv=0;
    end
end

%mDotThruOrifice calculates mass flow (mDot) in kg/s through an orifice
%given pressures and working fluid properties
function [mDot] = mDotThruOrifice(in1,in2,in3,in4,in5,in6)
%   Variable description:
%       in1: pressure upstream (Pa)
%       in2: pressure downstream (Pa)
%       in3: fluid density (kg/m^3)
%       in4: gamma = (fluidCP/fluidCV)
%       in5: orifice discharge coefficient
%       in6: orifice diameter (m)
%
%   Note on in1 and in2 variables:
%       The sign convention used in this function assumes normal (positive)
%       fluid movement from the in1 region to the in2 region (in1 is
%       upstream, in2 is downstream). However, these variables may be 
%       reversed so that there is flow from in2 to in1, however, the
%       resulting mDot will be negative.
%   
    if in1 < in2
        downP = in1;
        upP = in2;
        directionSign = -1;
    else
        upP = in1;
        downP = in2;
        directionSign = 1;
    end
    rho = in3;
    gamma = in4;
    outletCD = in5;
    outletDia = in6;
    outletArea = pi*(outletDia/2)^2;
    %fprintf("%0.3f, %0.3f, %0.6f, ", upP, downP, rho);
    %
    criticalP = upP * (2/(gamma+1))^(gamma/(gamma-1));
    %fprintf("%0.3f, %0.3f, ", criticalP, downP);
    if(downP < criticalP)
        %Choked
        %fprintf("C\n");
        r = downP/criticalP;
        r = (2/(gamma+1))^(gamma/(gamma-1));
        mDot = outletCD*outletArea*sqrt(upP*rho*(2*gamma/(gamma-1))*r^(2/gamma)*(1-r^((gamma-1)/gamma))); %kg/s
        %fprintf("%0.20f, ", rho);
    else
        %Subsonic
        %fprintf("S\n");
        r = downP/upP;
        mDot = outletCD*outletArea*sqrt(upP*rho*(2*gamma/(gamma-1))*r^(2/gamma)*(1-r^((gamma-1)/gamma))); %kg/s
    end
    %fprintf("\n");
    mDot = mDot * directionSign; %Corrects sign on mDot to follow stated convention above
end

%Outputs temperature [K], Pressure [Pa], and density [kg/m^3] given an altitude  
function [T, p, rho] = StandardAtm(h)
    
    rho_sea = 1.2250;     % density at sea level, kg/ m^3
    p_sea   = 1.01325E5;  % pressure at sea level, N/m^2
    R       = 287;        % gas constant, J/kg/K
    g_zero  = 9.81;       % gravitatoinal acceleration, m/s^2

    T_set   = [288.16,216.66,216.66,282.66,282.66,165.66,165.66,225.66]; % list of temperature points that define endpoints of each layer (starting at the ground), K
    h_set   = [0,11000,25000,47000,53000,79000,90000,105000]; %list of altitude points that define endpoints of each layer (starting at the ground), m
    a_set   = [-6.5*10^-3,0,3*10^-3,0,-4.5*10^-3,0,4*10^-3];%list of gradient layer slopes (starting at the ground), K/m
    p_set   = [p_sea, 2.27e4, 2.5273e3, 1.2558e2, 61.493, 1.2595, 0.162723]; % list of pressure at each layer endpoint, N/m^2
    rho_set = [rho_sea, 3.648e-1, 4.0639e-2, 1.5535e-3, 7.5791e-4, 2.184e-5, 1.84114e-6]; % list of density at each layer endpoint, kg/m^3
            
    if h <= h_set(2)
        layer = 1;
    elseif h <= h_set(3)
        layer = 2;
    elseif h <= h_set(4)
        layer = 3;
    elseif h <= h_set(5)
        layer = 4;
    elseif h <= h_set(6)
        layer = 5;
    elseif h <= h_set(7)
        layer = 6;
    %elseif h <= h_set(8)
    else
        layer = 7;
    end
    
    
    if mod(layer, 2) == 1
        %Gradient layer
        T = T_set(layer) + (a_set(layer) * (h - h_set(layer))); % temperature equation for gradient layer, K
        p = p_set(layer) * ((T / T_set(layer)) ^ (-g_zero / (a_set(layer) * R))); % pressure equation for gradient layer, N/m^2
        rho = rho_set(layer)* ((T / T_set(layer))^((-g_zero / (a_set(layer)*R))+1)); % density equation for gradient layer, kg/m^3
        
    else
        %Isothermal layers
        T = T_set(layer); % temperature for isothermal layer, K
        p = p_set(layer) * exp((-g_zero * (h - h_set(layer))) / (R * T)); % pressure equation for isothermal layer, N/m^2
        rho = (p * rho_set(layer)) / p_set(layer); % density equation for isothermal layer, kg/m^3
    end
    
end
       