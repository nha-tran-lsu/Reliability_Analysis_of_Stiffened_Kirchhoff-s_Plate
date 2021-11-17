function [Kp]=Assembling_Stiffened(option1,Emodule,Kp)
if option1 == 1
    [ Ksx ] = cal_Ksx(Emodule);
    [ Tsx, Tsy ] = Assembling_Ksx_Ksy;
    Kp=Kp+Tsx'*Ksx*Tsx;
elseif option1 == 2
    [ Ksx ] = cal_Ksx(Emodule);
    [ Ksy ] = cal_Ksy(Emodule);
    [ Tsx, Tsy ] = Assembling_Ksx_Ksy;
    Kp = Kp + Tsx'*Ksx*Tsx + Tsy'*Ksy*Tsy ;
end