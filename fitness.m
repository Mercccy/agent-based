function y = fitness(z)
    global x equation
    x=z;
    y =abs(eval(equation));
end