function [ y ] = compute_tanh( x )
y=(exp(x)-exp(-x))/(exp(x)+exp(-x));
end

