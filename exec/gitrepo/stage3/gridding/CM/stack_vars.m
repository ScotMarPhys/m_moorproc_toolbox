function Ustacked = stack_vars(U1,U2,U3,U4,U5,U6,U7)
    Ustacked = cat(3,U1,U2,U3,U4,U5,U6,U7);
    Ustacked = sum(Ustacked,3,'omitnan');
    Ustacked(Ustacked==0)=NaN;
end