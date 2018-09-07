function y=trascendent_term(x)
indMin=find(x<=731);
indMax =find(x>731);
y(indMin)=((pi*x(indMin)).^(-1/2).*exp(-x(indMin))./(erfc(x(indMin).^(1/2))));   
y(indMax)=1;    
end
