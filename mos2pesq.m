function p=mos2pesq(m)
%mos is MOS-LQO calculated by PESQ C code
a=0.999;
b=4.999-a;
c=-1.4945;
d=4.6607;
p=(log(b./(m-a)-1)-d)/c;
end
