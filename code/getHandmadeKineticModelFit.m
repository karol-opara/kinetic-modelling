function kHandmade = getHandmadeKineticModelFit()
k4=[0.01
    0.10
    0.03
    0.01
    0.06
    0.03];

k5=[0.01
    0.03
    0.08
    0.01
    0.3
    0.07];

k7=[0.025
    0.06
    0.08
    0.03
    0.11
    0.025];

k6=[0.027
    0.08
    0.07
    0.03
    0.13
    0.03];
   
kHandmade = {};
kHandmade{4}=k4;
kHandmade{5} =k5;
kHandmade{6}=k6;
kHandmade{7}=k7;
end