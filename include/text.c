if(k==0&& l == 0) c=((power(h2,2) + power(h3,2))*(3 + cos(2*Jpsithe)))/2. + 2*power(h1,2)*power(sin(Jpsithe),2);

if(k==0&& l == 1) c=-(sqrt(0.4)*h1*sin(2*Jpsithe)*(2*sqrt(3)*h2*sin(phi2) + 3*h3*sin(phi3)))/3.;

if(k==0&& l == 6) c=(-((h2 - h3)*(h2 + h3)*(3 + cos(2*Jpsithe))) - 4*power(h1,2)*power(sin(Jpsithe),2))/(2.*sqrt(3));

if(k==0&& l == 7) c=sqrt(2.0/3)*h1*h3*cos(phi3)*sin(2*Jpsithe);

if(k==0&& l == 8) c=(2*h2*h3*cos(phi2 - phi3)*power(sin(Jpsithe),2))/sqrt(3);

if(k==0&& l == 10) c=(2*h2*h3*power(sin(Jpsithe),2)*sin(phi2 - phi3))/sqrt(3);

if(k==0&& l == 11) c=(2*h1*sin(2*Jpsithe)*(3*h2*sin(phi2) - sqrt(3)*h3*sin(phi3)))/(3.*sqrt(5));

if(k==1&& l == 2) c=-(sqrt(2.0/15)*h1*h2*cos(phi2)*sin(2*Jpsithe));

if(k==1&& l == 3) c=-((h2*h3*(3 + cos(2*Jpsithe))*cos(phi2 - phi3))/sqrt(5)) - (2*(2*power(h1,2) + power(h2,2))*power(sin(Jpsithe),2))/sqrt(15);

if(k==1&& l == 4) c=sqrt(2.0/3)*h1*h3*sin(2*Jpsithe)*sin(phi3);

if(k==1&& l == 5) c=-((h2*h3*(3 + cos(2*Jpsithe))*sin(phi2 - phi3))/sqrt(3));

if(k==1&& l == 12) c=sqrt(1.2)*h1*h2*cos(phi2)*sin(2*Jpsithe);

if(k==1&& l == 13) c=-(sqrt(2.0/15)*h2*h3*(3 + cos(2*Jpsithe))*cos(phi2 - phi3)) + sqrt(0.4)*(2*power(h1,2) + power(h2,2))*power(sin(Jpsithe),2);

if(k==1&& l == 14) c=-(sqrt(2.0/3)*h1*h3*cos(phi3)*sin(2*Jpsithe));

if(k==1&& l == 15) c=-(sqrt(2.0/3)*power(h3,2)*power(sin(Jpsithe),2));

if(k==2&& l == 0) c=-2*sqrt(2)*h1*h2*cos(Jpsithe)*sin(Jpsithe)*sin(phi2);

if(k==2&& l == 1) c=(h2*h3*(3 + cos(2*Jpsithe))*cos(phi2 - phi3))/sqrt(5) - (2*(-2*power(h1,2) + power(h2,2))*power(sin(Jpsithe),2))/sqrt(15);

if(k==2&& l == 6) c=sqrt(2.0/3)*h1*h2*sin(2*Jpsithe)*sin(phi2);

if(k==2&& l == 7) c=-((h2*h3*(3 + cos(2*Jpsithe))*sin(phi2 - phi3))/sqrt(3));

if(k==2&& l == 8) c=sqrt(2.0/3)*h1*h3*sin(2*Jpsithe)*sin(phi3);

if(k==2&& l == 9) c=sqrt(2.0/3)*power(h3,2)*power(sin(Jpsithe),2);

if(k==2&& l == 10) c=sqrt(2.0/3)*h1*h3*cos(phi3)*sin(2*Jpsithe);

if(k==2&& l == 11) c=sqrt(2.0/15)*h2*h3*(3 + cos(2*Jpsithe))*cos(phi2 - phi3) + sqrt(0.4)*(-2*power(h1,2) + power(h2,2))*power(sin(Jpsithe),2);

if(k==3&& l == 2) c=(-((power(h2,2) - 3*power(h3,2))*(3 + cos(2*Jpsithe))) + 4*power(h1,2)*power(sin(Jpsithe),2))/(2.*sqrt(15));

if(k==3&& l == 3) c=(sqrt(0.4)*h1*(-2*sqrt(3)*h2*cos(phi2) + 3*h3*cos(phi3))*sin(2*Jpsithe))/3.;

if(k==3&& l == 4) c=(2*h2*h3*power(sin(Jpsithe),2)*sin(phi2 - phi3))/sqrt(3);

if(k==3&& l == 5) c=-(sqrt(2.0/3)*h1*h3*sin(2*Jpsithe)*sin(phi3));

if(k==3&& l == 12) c=((3*power(h2,2) + power(h3,2))*(3 + cos(2*Jpsithe)) - 12*power(h1,2)*power(sin(Jpsithe),2))/(2.*sqrt(15));

if(k==3&& l == 13) c=(2*h1*(3*h2*cos(phi2) + sqrt(3)*h3*cos(phi3))*sin(2*Jpsithe))/(3.*sqrt(5));

if(k==3&& l == 14) c=(2*h2*h3*cos(phi2 - phi3)*power(sin(Jpsithe),2))/sqrt(3);


















if(k==0&& l == 0) c=1;

if(k==0&& l == 3) c=alpha;

if(k==1&& l == 0) c=alpha*cos(phi)*sin(theta);

if(k==1&& l == 1) c=gamma*cos(phi)*cos(theta) - beta*sin(phi);

if(k==1&& l == 2) c=-(beta*cos(phi)*cos(theta)) - gamma*sin(phi);

if(k==1&& l == 3) c=cos(phi)*sin(theta);

if(k==2&& l == 0) c=alpha*sin(phi)*sin(theta);

if(k==2&& l == 1) c=beta*cos(phi) + gamma*cos(theta)*sin(phi);

if(k==2&& l == 2) c=gamma*cos(phi) - beta*cos(theta)*sin(phi);

if(k==2&& l == 3) c=sin(phi)*sin(theta);

if(k==3&& l == 0) c=alpha*cos(theta);

if(k==3&& l == 1) c=-(gamma*sin(theta));

if(k==3&& l == 2) c=beta*sin(theta);

if(k==3&& l == 3) c=cos(theta);




if(k==0&& l == 0) c=1;

if(k==0&& l == 3) c=alpha;

if(k==1&& l == 0) c=sqrt(0.6)*alpha*sin(phi)*sin(theta);

if(k==1&& l == 1) c=2*sqrt(0.6)*(beta*cos(phi) + gamma*cos(theta)*sin(phi));

if(k==1&& l == 2) c=2*sqrt(0.6)*(gamma*cos(phi) + beta*cos(theta)*sin(phi));

if(k==1&& l == 3) c=sqrt(0.6)*sin(phi)*sin(theta);

if(k==6&& l == 0) c=-(sqrt(3)*(1 + 3*cos(2*theta)))/4.;

if(k==6&& l == 3) c=-(sqrt(3)*alpha*(1 + 3*cos(2*theta)))/4.;

if(k==7&& l == 0) c=-3*cos(phi)*cos(theta)*sin(theta);

if(k==7&& l == 3) c=-3*alpha*cos(phi)*cos(theta)*sin(theta);

if(k==8&& l == 0) c=(-3*cos(2*phi)*Power(sin(theta),2))/2.;

if(k==8&& l == 3) c=(-3*alpha*cos(2*phi)*Power(sin(theta),2))/2.;

if(k==10&& l == 0) c=-9*alpha*cos(phi)*cos(theta)*sin(phi)*Power(sin(theta),2);

if(k==10&& l == 1) c=-3*beta*cos(2*phi)*cos(theta)*sin(theta) - (3*gamma*(1 + 3*cos(2*theta))*sin(2*phi)*sin(theta))/4.;

if(k==10&& l == 2) c=-3*gamma*cos(2*phi)*cos(theta)*sin(theta) + (3*beta*(1 + 3*cos(2*theta))*sin(2*phi)*sin(theta))/4.;

if(k==10&& l == 3) c=-9*cos(phi)*cos(theta)*sin(phi)*Power(sin(theta),2);

if(k==11&& l == 0) c=(-9*alpha*(3 + 5*cos(2*theta))*sin(phi)*sin(theta))/(4.*sqrt(10));

if(k==11&& l == 1) c=(-3*beta*cos(phi)*(3 + 5*cos(2*theta)))/(4.*sqrt(10)) - (3*gamma*(cos(theta) + 15*cos(3*theta))*sin(phi))/(8.*sqrt(10));

if(k==11&& l == 2) c=(-3*gamma*cos(phi)*(3 + 5*cos(2*theta)))/(4.*sqrt(10)) + (3*beta*(cos(theta) + 15*cos(3*theta))*sin(phi))/(8.*sqrt(10));

if(k==11&& l == 3) c=(-9*(3 + 5*cos(2*theta))*sin(phi)*sin(theta))/(4.*sqrt(10));
















