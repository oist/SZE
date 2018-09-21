function[symmetry] = checkSymmetry2D(state2D)
    normalDev = max(max(state2D - transpose(state2D)))/max(state2D(:));
    symmetry = normalDev<1;
end