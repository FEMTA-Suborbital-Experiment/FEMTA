clear
clc
close all
%https://multimedia.3m.com/mws/media/199819O/3m-novec-7200-engineered-fluid-en.pdf

%"tank" suffix denotes propellant tank values
%"CC" suffix denotes collection chamber values
%"O" suffix denotes orifice values
%"pipe" suffix denotes piping network values

%CONSTANT AND INITIAL VALUES
global R kB N_a;
R=8.3145;           %Universal Gas Constant [kJ/kmol-K]
kB=1.380649e-23;    %Boltzmann constant [J/K]
N_a=6.02214e23;     %Avogadros number [# of particles/mol]

%Air Properties
gammaAir=1.4;     %ratio of specific heats for air
MW_Air=29;        %Molecular Weight of Air [kg/kmol]
R_air=R/MW_Air;   %Specific Gas Constant for Air [kJ/kg-K]
Cp_Air=1.005;     %Specific Heat of Air [kJ/kg-K]
 
%Water Properties
rho_water=997;      %density of liquid water [kg/m^3]
m_H2O=2.988e-26;    %mass of one water molecule [kg]
MW_Water=18.0135;   %Molecular weight of water [kg/kmol]
h_evap_water=2.260; %heat of vaporazitaion of water [kJ/kg]
Cp_water_liquid=4;  %specific heat of liquid water [kJ/kg-K]
Cp_water_vapor=2;   %specific heat of water vapor [kJ/kg-K]
Ce_water=1;         %Evaporation Coefficient of water
Cc_water=0.75;      %Condensation Coefficient of water

%HFE Properties
MW_HFE=250;               %Molecular Weight of HFE [kg/kmol]
R_HFE=R/MW_HFE;           %Specific Gas Constant for HFE [kJ/kg-K]
h_evap_HFE=0.0308/MW_HFE; %Heat of Vaporization of HFE [kJ/kg]
Cp_HFEliquid=1.172303;    %Specific Heat of HFE liquid [kJ/kg-K]
m_HFE=MW_HFE/(N_a*1000);  %mass of one HFE molecule [kg]
Ce_HFE=1;          %Evaporation Coefficient of HFE
Cc_HFE=0.1;         %Condensation Coefficient of HFE

%Propellant Tank
P0_tank=101325;                                     %Initial pressure in prop tank [Pa] (1 atm)
T0_tank=300;                                        %Initial temperature in prop tank [K]
V_tank=3e-5;                                        %Total volume of prop tank [m^3] (30mL) (CONSTANT)
volWater0_tank=1.4e-5;                              %initial volume of water in prop tank [m^3] (14mL)
volHFE_liquid0_tank=15e-6;                          %initial volume of HFE in prop tank [m^3] (10mL)
volAir0=V_tank-volWater0_tank-volHFE_liquid0_tank;  %Initial volume of air in prop tank [m^3]
n_Air=(P0_tank*volAir0)/(R*T0_tank);                %Number of moles of air in propr tank (CONSTANT)
A_HFE=3.167e-5;                                     %Area from which HFE condenses and evaporates [m^2]

%Collection Chamber
P0_CC=0;        %initial pressure in collection chamber [Pa]
T0_CC=300;      %initial temperature in collection chamber [K]
V_CC=227.75e-6; %volume of CC [m^3] (227.75mL)

%Piping network
D_pipe=(1/8)/39.37;  %Pipe diameter [m] (CONSTANT)
A=pi*((D_pipe/2)^2); %Area of pipe cross section [m^2] (CONSTANT)

%Orifice
D_O=0.000127;          %Orifice diameter [m]
A_O=pi*((D_O/2).^2); %Area of orifice [m^2]
Beta=D_O/D_pipe;     %Ratio of orifice to pipe diameter
CD_orifice=0.6;      %discharge coefficient of orifice

%VARIABLES
volWater_tank=volWater0_tank;      %Initial volume of water in Prop Tank [m^3]
volWater_CC=0;                     %Initial volume of water in Collection Chamber [m^3]
volHFE_liquid=volHFE_liquid0_tank; %Volume of HFE Liquid in Prop Tank [m^3]
tankPress=P0_tank;                 %Pressure in Prop Tank [Pa]
CCPress=P0_CC;                    %Pressure in Collection Chamber [Pa]
tankTempGas=T0_tank;               %Temperature of gas in Prop Tank [K]
tankTempLiquid=T0_tank;            %Temperature of HFE Liquid in Prop Tank [K]
CCTempGas=T0_CC;                   %Temperature of gas in CC [K]
CCTempLiquid=T0_CC;                %Temperature of Water Liquid in CC [K]
m_HFE_vapor=0;                     %Initial mass of HFE vapor in prop tank
m_HFE_liquid=volHFE_liquid0_tank*nvcRho(T0_tank);  %Initial mass of HFE liquid in prop tank
m_water_vapor=0;                   %Initial mass of water vapor in CC
m_water_liquid=0;                  %Initial mass of water liquid in CC
n_Gas=n_Air;
time=0;
count=1;
dt=1e-3;                           %timestep [s]

%ARRAY INITIALIZATION
tankVolWater_array=[];
tankPress_array=[];
tankTempGas_array=[];
tankTempLiquid_array=[];
CCVolWater_array=[];
CCPress_array=[];
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

%EXECUTABLE LOOP
while volWater_tank>0
    
    %Volumetric Flow Rate of Liquid Water Propellant [m^3/s] 
    flo_water=CD_orifice*A_O*sqrt(2*(tankPress-CCPress)/(rho_water*(1-(Beta^4))));
    volWater_tank=volWater_tank-(flo_water*dt);
    volWater_CC=volWater_CC+(2*flo_water*dt);
    
    %PROPELLANT TANK
        %Vapor Pressure of HFE [Pa]
        Pvap_HFE=nvcVP(tankTempLiquid);

        %Density of HFE liquid [kg/m^3]
        rho_HFE=nvcRho(tankTempLiquid);
        
        %Volume, Temperature, and Pressure of Gas
        volGas=V_tank-volWater_tank-volHFE_liquid;
        tankPress=(n_Gas*R*tankTempGas)/volGas;
        %tankTempGas=T0_tank*((volAir0/volGas)^(gammaAir-1));
        
        %Amount of HFE either condensing or evaporating
        m_HFE_transfer=HerKnu(Pvap_HFE,tankTempLiquid,tankTempGas,tankPress,m_HFE,A_HFE,Ce_HFE,Cc_HFE)*dt;
        
        %No mass evaporating if no liquid and no mass condensing if no vapor
        if (m_HFE_vapor==0 && m_HFE_transfer<0) || (m_HFE_liquid==0 && m_HFE_transfer>0)
            m_HFE_transfer=0;
        end
        
        %Mass of vapor and Liquid is zero if negative
        if m_HFE_liquid<0
            m_HFE_liquid=0;
        elseif m_HFE_vapor<0
            m_HFE_vapor=0;
        end
        
        if tankPress>Pvap_HFE && m_HFE_transfer>0
            m_HFE_transfer=0;
        elseif tankPress>Pvap_HFE && m_HFE_transfer<0
            if m_HFE_vapor==0
                m_HFE_transfer=0;
            end
        end
        
        %Update amount of HFE liquid 
        m_HFE_liquid=m_HFE_liquid-m_HFE_transfer;
        vol_HFE_liquid=m_HFE_liquid/rho_HFE;
        
        %Update amount of HFE vapor 
        m_HFE_vapor=m_HFE_transfer+m_HFE_vapor;
        n_HFE_vapor=m_HFE_vapor/MW_HFE;
        
        %Update total amount of Gas (Air + HFE)
        n_Gas=n_HFE_vapor+n_Air;
        volGas=V_tank-volWater_tank-vol_HFE_liquid;
        
        %Temperatrue Update
        if m_HFE_liquid>0
            Q_HFE=m_HFE_transfer*h_evap_HFE;
            tankTempLiquid=(-Q_HFE/(m_HFE_liquid*Cp_HFEliquid))+tankTempLiquid;
        end
        
    %COLLECTION CHAMBER
        %Surface Area of Collected Water [m^2]
        r=((3*volWater_CC)/(4*pi))^(1/3);
        A_water=4*pi*(r^2);
        
        %Vapor Pressure of Water [Pa]
        Pvap_water=waterVP(CCTempLiquid);
        
        %Mass of water either evaporating or condensing at current timestep [kg]
        m_water_transfer=HerKnu(Pvap_water,CCTempLiquid,CCTempGas,CCPress,m_H2O,A_water,Ce_water,Cc_water)*dt;
        
        %No mass evaporating if no liquid and no mass condensing if no vapor
        if (m_water_vapor==0 && m_water_transfer<0) || (m_water_liquid==0 && m_water_transfer>0)
            m_HFE_transfer=0;
        end
        
        if CCPress>Pvap_water && m_water_transfer>0
            m_water_transfer=0;
        elseif CCPress>Pvap_HFE && m_water_transfer<0
            if m_water_vapor==0
                m_water_transfer=0;
            end
        end
        
        %Total mass of water vapor at current time [kg]
        m_water_vapor=m_water_transfer+m_water_vapor;
        
        %Total mass of water liquid at current time [kg]
        m_water_liquid=(volWater_CC*rho_water)-m_water_vapor;
        
        %Mass of Vapor and Liquid is zero if negative
        if m_water_liquid<0
            m_water_liquid=0;
        elseif m_water_vapor<0
            m_water_vapor=0;
        end
        
        n_water_vapor=m_water_vapor/MW_Water;
        
        %Pressure Update
        CCPress=(n_water_vapor*R*CCTempGas)/(V_CC-volWater_CC);
        
        %Temperature Update
        Q_water=m_water_transfer*h_evap_water;
        
        if m_water_liquid>0
            CCTempLiquid=(-Q_water/(m_water_liquid*Cp_water_liquid))+CCTempLiquid;
        end
    
    %ARRAY UPDATE
        %Propellant Tank
        tankVolWater_array(count)=volWater_tank;
        tankPress_array(count)=tankPress;
        tankTempGas_array(count)=tankTempGas;
        tankTempLiquid_array(count)=tankTempLiquid;
        tankVolGas_array(count)=volGas;
        m_HFE_vapor_array(count)=m_HFE_vapor;
        m_HFE_liquid_array(count)=m_HFE_liquid;
        n_Gas_array(count)=n_Gas;
        
        %Collection Chamber
        CCVolWater_array(count)=volWater_CC;
        CCPress_array(count)=CCPress;
        CCTempGas_array(count)=CCTempGas;
        CCTempLiquid_array(count)=CCTempLiquid;
        m_water_vapor_array(count)=m_water_vapor;
        m_water_liquid_array(count)=m_water_liquid;
        m_water_transfer_array(count)=m_water_transfer;
        
        %Miscellaneous
        PvapHFE_array(count)=Pvap_HFE;
        QHFE_array(count)=Q_HFE;
        PvapWater_array(count)=Pvap_water;
        Qwater_array(count)=flo_water;
        time_array(count)=time;
        m_HFE_transfer_array(count)=m_HFE_transfer;
        
        time=time+dt;
        count=count+1;      
end

%Data manipulation
time_array=time_array';

%Plots
figure(1)
plot(time_array,tankVolWater_array'*10^6,'Linewidth',3)
hold on
plot(time_array,CCVolWater_array'*10^6,'Linewidth',3)
xlabel("Time [s]",'Fontsize',17)
ylabel("Volume of Water in Prop Tank and CC [mL]",'Fontsize',17)
legend("Tank","CC",'Fontsize',15)
title("Volumes of Collection Chamber and Propellant Tank",'Fontsize',22)

figure(2)
plot(time_array,tankTempGas_array')
hold on
plot(time_array,CCTempGas_array')
plot(time_array,CCTempLiquid_array')
plot(time_array,tankTempLiquid_array')
xlabel("time [s]",'Fontsize',17)
ylabel("Temperature [K]",'Fontsize',17)
legend("Temp of gas in tank","Temp of gas in CC","temp of liquid in CC","Temp of liquid in tank",'Fontsize',17)
title("Temperature values for Duration of Experiment",'Fontsize',22)

figure(3)
plot(time_array,tankPress_array'/1000,'Linewidth',3)
hold on
plot(time_array,PvapHFE_array'/1000,'Linewidth',3)
xlabel("Time [s]",'Fontsize',17)
ylabel("Pressure of Prop Tank [kPa]",'Fontsize',17)
legend("Tank Pressure","Vapor Pressure of HFE",'Fontsize',17)
title("Pressures in Propellant Tank",'Fontsize',22)

figure(4)
plot(time_array,CCPress_array'/1000,'Linewidth',3)
hold on
plot(time_array,PvapWater_array'/1000,'Linewidth',3)
xlabel("Time [s]",'Fontsize',17)
ylabel("Pressure of Collection Chamber [kPa]",'Fontsize',17)
legend("CC Pressure","Vapor Pressure of Water",'Fontsize',17)
title("Pressures in Collection Chamber",'Fontsize',22)

figure(5)
plot(time_array,m_water_transfer_array')


% figure(3)
% plot(timescale,tankTemp_array)
% xlabel("time [s]")
% ylabel("Temperature of Prop Tank [K]")

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
function m_transfer=HerKnu(Ps,T_liquid,T_vapor,Pg,m,A,C_evap,C_cond)
    global kB;
    m_transfer=A*(sqrt(m/(2*pi*kB))*((C_evap*(Ps/sqrt(T_liquid)))-(C_cond*(Pg/sqrt(T_vapor)))));
end


