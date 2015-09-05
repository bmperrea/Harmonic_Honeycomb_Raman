//author: Brent Perreault
//credit: based off of an images by Philippe Ivald, Charles Staats, as well as one by g.kov, all on tex.stackexchange.com
//see questions 141363 and 121169 for originals. See also the notes by Staats: //http://math.uchicago.edu/~cstaats/Charles_Staats_III/Notes_and_papers_files/asymptote_tutorial.pdf

import three;
settings.render=8;
settings.prc=false;
size(10cm);

//origin
O=(0,0,0);

//perspective
triple persp = (-6,1.5,-1);
real n = 4;
currentprojection = perspective(n*persp,up=Y);

//Initial data

pen fontp = fontsize(14pt);
pen unitp = black+1.5bp;
pen framep=gray(0.4)+dashed+1.3bp;

material sphereCcolor = material(ambientpen= yellow, diffusepen=yellow, emissivepen=gray(.3));
material veccolor = material(diffusepen=red, ambientpen=red);
material framecolor = material(diffusepen=white,  ambientpen=white);
material sphereCcolor2 = material(ambientpen=gray(.5), diffusepen=gray(.5),  emissivepen=gray(.3));

material cylcolorz = material(diffusepen=blue, ambientpen=deepblue);
material cylcolory = material(diffusepen=green, ambientpen=deepgreen);
material cylcolorx = material(diffusepen=red, ambientpen=deepred);
material cylcolorg = material(diffusepen=lightgray, ambientpen=darkgray);

real cylRadius = 0.033;
real sphereRadius = 0.2;
real frameCyl = .015;
real frameRadius = frameCyl;

//In these functions the smaller z value should come first

void drawRod(triple a, triple b, material cylcolor) {
  surface rod = extrude(scale(cylRadius)*unitcircle, axis=length(b-a)*Z);
  triple orthovector = cross(Z, b-a);
  if (length(orthovector) > .01) {
    real angle = aCos(dot(Z, b-a) / length(b-a));
    rod = rotate(angle, orthovector) * rod;
  }
  draw(shift(a)*rod, surfacepen=cylcolor);
}

void drawFrame(triple a, triple b) {
	surface rod = extrude(scale(frameCyl)*unitcircle, axis=length(b-a)*Z);
	triple orthovector = cross(Z, b-a);
	if (length(orthovector) > .01) {
		real angle = aCos(dot(Z, b-a) / length(b-a));
		rod = rotate(angle, orthovector) * rod;
	}
	draw(shift(a)*rod, surfacepen=framecolor);
	draw(shift(b)*scale3(frameRadius)*unitsphere, surfacepen=framecolor);
}

void drawCarbon(triple center) {
draw(shift(center)*scale3(sphereRadius)*unitsphere, surfacepen=sphereCcolor);
}

void drawCarbon2(triple center) {
     draw(shift(center)*scale3(sphereRadius)*unitsphere, surfacepen=sphereCcolor2);
}

// lattice vectors from Mandal and Surendran PRB 79 024426

triple dz = (0,0,sqrt(2));
triple dxa = (-1/sqrt(3),-sqrt(2/3),-1);
triple dxb = (-1/sqrt(3),sqrt(2/3),-1);
triple dya = (1/sqrt(3),sqrt(2/3),-1);
triple dyb = (1/sqrt(3),-sqrt(2/3),-1);

//Lattice sites

//The sign for dx depends on sublattice which switches with every hop
//The sign for dz is always negative going from Aa and always postive out the other way (Ab)
//The choice dxa or dxb (same for y) depends on which ladder you are on, which switches with every z-hop

triple Aa = O; 
triple Ab = dxa;  
//draw the yellow ones first because the first thing appears to survive.
drawCarbon(Aa); //consider this on the 2 site (b sublattice)
drawCarbon(Ab);
drawRod(Aa,Ab, material(diffusepen=yellow, ambientpen=lightolive)); 

transform3 tran = rotate(120,dxa,dz);
triple move = tran*O;
tran = shift(dxa/2-move)*tran;

//An arc around the x-bond
draw(tran*arc(O,.5,90,-50,90,300,(0,0,1)),black+linewidth(1pt),ArcArrow3);

//oneflux
triple Aaz = Aa + dz; drawCarbon2(Aaz); drawRod(Aa,Aaz,cylcolorz); 
triple Aazy = Aaz - dyb; drawCarbon2(Aazy); drawRod(Aazy,Aaz,cylcolory); 
triple Aazyx = Aazy + dxb; drawCarbon2(Aazyx); drawRod(Aazyx,Aazy,cylcolorx); 
triple Aazyxz = Aazyx - dz; drawCarbon2(Aazyxz); drawRod(Aazyxz,Aazyx,cylcolorz);

triple Abz = Ab - dz; drawCarbon2(Abz); drawRod(Abz,Ab,cylcolorz); 
triple Abzx = Abz + dxb; drawCarbon2(Abzx); drawRod(Abzx,Abz,cylcolorx); 
triple Abzxy = Abzx - dyb; drawCarbon2(Abzxy); drawRod(Abzxy,Abzx,cylcolory); 
triple Abzxyz = Abzxy + dz; drawCarbon2(Abzxyz); drawRod(Abzxy,Abzxyz,cylcolorz);

drawRod(Aa,Ab,cylcolorx);
drawRod(Abzxyz,Aazyxz,cylcolorx);
drawCarbon2(Aa); //consider this on the 2 site (b sublattice)
drawCarbon2(Ab);

//translate
triple t = Aa-Aazyxz;

triple Aa2 = Aa+t;
triple Aaz2 = Aaz+t;
triple Aazy2 = Aazy + t;
triple Aazyx2 = Aazyx + t;
triple Aazyxz2 = Aazyxz + t;
triple Ab2 = Ab+t;
triple Abz2 = Abz+t;
triple Abzx2 = Abzx + t;
triple Abzxy2 = Abzxy + t;
triple Abzxyz2 = Abzxyz + t;

//twoflux

drawCarbon2(Aaz2); drawRod(Aa2,Aaz2,cylcolorz);
drawCarbon2(Aazy2); drawRod(Aazy2,Aaz2,cylcolory);
drawCarbon2(Aazyx2); //drawRod(Aazyx2,Aazy2,cylcolorx); 
drawCarbon2(Aazyxz2); drawRod(Aazyxz2,Aazyx2,cylcolorz);
drawCarbon2(Abz2); drawRod(Abz2,Ab2,cylcolorz); 
drawCarbon2(Abzx2); drawRod(Abzx2,Abz2,cylcolorx);
drawCarbon2(Abzxy2); drawRod(Abzxy2,Abzx2,cylcolory);
drawCarbon2(Abzxyz2); drawRod(Abzxy2,Abzxyz2,cylcolorz);

drawRod(Aa2,Ab2,cylcolorx);
//drawRod(Abzxyz2,Aazyxz2,cylcolorx);
drawCarbon2(Aa2); //consider this on the 2 site (b sublattice)
drawCarbon2(Ab2);



//Three flux (from scratch)
triple Aay = Ab - dya; drawCarbon2(Aay); drawRod(Aay,Ab,cylcolory);
triple Aayz = Aay + dz; drawCarbon2(Aayz); drawRod(Aay,Aayz,cylcolorz);
triple Aayzx = Aayz - dxb; drawCarbon2(Aayzx); drawRod(Aayzx,Aayz,cylcolorx);
triple Aayzxz = Aayzx + dz; drawCarbon2(Aayzxz); drawRod(Aayzx,Aayzxz,cylcolorz);

triple Abz = Aa + dz; drawCarbon2(Abz); drawRod(Aa,Abz,cylcolorz);
triple Abzx = Abz - dxb; drawCarbon2(Abzx); drawRod(Abz,Abzx,cylcolorx);
triple Abzxz = Abzx + dz; drawCarbon2(Abzxz); drawRod(Abzx,Abzxz,cylcolorz);
triple Abzxzy = Abzxz - dya; drawCarbon2(Abzxzy); drawRod(Abzxzy,Abzxz,cylcolory);

drawRod(Abzxzy,Aayzxz,cylcolorx);


//Four flux (from scratch)
triple Aay = Ab - dya; drawCarbon2(Aay); drawRod(Aay,Ab,cylcolory);
triple Aayz = Aay + dz; drawCarbon2(Aayz); drawRod(Aay,Aayz,cylcolorz);
triple Aayzy = Aayz - dyb; drawCarbon2(Aayzy); drawRod(Aayzy,Aayz,cylcolory);
triple Aayzyz = Aayzy + dz; drawCarbon2(Aayzyz); drawRod(Aayzy,Aayzyz,cylcolorz);

triple Abz = Aa + dz; drawCarbon2(Abz); drawRod(Aa,Abz,cylcolorz);
triple Abzy = Abz - dyb; drawCarbon2(Abzy); drawRod(Abz,Abzy,cylcolory);
triple Abzyz = Abzy + dz; drawCarbon2(Abzyz); drawRod(Abzy,Abzyz,cylcolorz);
triple Abzyzy = Abzyz - dya; drawCarbon2(Abzyzy); drawRod(Abzyzy,Abzyz,cylcolory);

drawRod(Abzyzy,Aayzyz,cylcolorx);

//Flux five (from scratch)
triple Aay = Aa + dya; drawCarbon2(Aay); drawRod(Aay,Aa,cylcolory); 
triple Aayz = Aay - dz; drawCarbon2(Aayz); drawRod(Aayz,Aay,cylcolorz); 
triple Aayzy = Aayz + dxb; drawCarbon2(Aayzy); drawRod(Aayzy,Aayz,cylcolorx); 
triple Aayzyz = Aayzy - dz; drawCarbon2(Aayzyz); drawRod(Aayzyz,Aayzy,cylcolorz);

triple Abz = Ab - dz; drawCarbon2(Abz); drawRod(Abz,Ab,cylcolorz); 
triple Abzy = Abz + dxb; drawCarbon2(Abzy); drawRod(Abz,Abzy,cylcolorx); 
triple Abzyz = Abzy - dz; drawCarbon2(Abzyz); drawRod(Abzyz,Abzy,cylcolorz); 
triple Abzyzy = Abzyz + dya; drawCarbon2(Abzyzy); drawRod(Abzyz,Abzyzy,cylcolory);

drawRod(Abzyzy,Aayzyz,cylcolorx);


//Flux six (from scratch)
triple Aay = Aa + dya; drawCarbon2(Aay); drawRod(Aay,Aa,cylcolory); 
triple Aayz = Aay - dz; drawCarbon2(Aayz); drawRod(Aayz,Aay,cylcolorz); 
triple Aayzx = Aayz + dyb; drawCarbon2(Aayzx); drawRod(Aayzx,Aayz,cylcolory); 
triple Aayzxz = Aayzx - dz; drawCarbon2(Aayzxz); drawRod(Aayzxz,Aayzx,cylcolorz);

triple Abz = Ab - dz; drawCarbon2(Abz); drawRod(Abz,Ab,cylcolorz); 
triple Abzx = Abz + dyb; drawCarbon2(Abzx); drawRod(Abz,Abzx,cylcolory); 
triple Abzxz = Abzx - dz; drawCarbon2(Abzxz); drawRod(Abzxz,Abzx,cylcolorz); 
triple Abzxzy = Abzxz + dya; drawCarbon2(Abzxzy); drawRod(Abzxz,Abzxzy,cylcolory);

drawRod(Abzxzy,Aayzxz,cylcolorx);




triple Ac = Ab - dxa;
triple Ad = Ac + dz;

triple Ba = Ad - dyb;
triple Ca = Ad - dxb; 
triple Db = Ac + dya;
triple Da = Db - dz;

triple Ea = (2 xpart(Da),0,0);
triple Eb = Ea + dz;
triple Fa = Ba + Ea;

//Unit Vectors
triple aa = Ca;
triple ab = Ba;
triple ac = Da;



//Alternative unit vectors that are useful
triple ad = -ab+aa;   //(goes along the ladder that a3 does not, doing so also with a positive b-component)
triple ae = ac - ad;  //(lies in a-direction)




draw(Aa--aa,arrowhead,themargin,L=Label("$a_2$",2aa+1.5Y-2Z,p=fontp),p=unitp);
draw(Aa--ab,arrowhead,themargin,L=Label("$a_1$",2ab+3Y+2.2Z,p=fontp),p=unitp);
draw(Aa--ac,arrowhead,themargin,L=Label("$a_3$",2ac-4Y-.6Z,p=fontp),p=unitp);

//dimensions of the unit cell
real uX = xpart(Ea);
real uY = ypart(Ba);
real uZ = zpart(Ba);





//abc axes

Label lx=Label("$a$",-1.2Z,p=fontp);
Label ly=Label("$b$",Z+.9Y,p=fontp);
Label lz=Label("$c$",.95Y+1.2Z,p=fontp);

triple axesOrigin = (-1.5,2,1.6);
arrowbar3 arrowaxes = Arrow3(size=7);

draw( shift(axesOrigin)*( O--X/2 ), L=lx ,arrow=arrowaxes,p=black+1bp); 
draw( shift(axesOrigin)*( O--Y/2 ) , L=ly, arrow=arrowaxes,p=black+1bp); 
draw( shift(axesOrigin)*( O--Z/2 ) , L=lz, arrow=arrowaxes,p=black+1bp); 


//xyz axes

Label lx=Label("$x$",.8Y+1.6Z,p=fontp);
Label ly=Label("$y$",3Z+.6Y,p=fontp);
Label lz=Label("$z$",1.2Y+Z,p=fontp);

triple axesOrigin = (-1.5,2,2.5);
arrowbar3 arrowaxes = Arrow3(size=7);

draw( shift(axesOrigin)*( O--(Z-X)/sqrt(8) ), L=lx ,arrow=arrowaxes,p=black+.5bp); 
draw( shift(axesOrigin)*( O--(X+Z)/sqrt(8) ) , L=ly, arrow=arrowaxes,p=black+.5bp); 
draw( shift(axesOrigin)*( O--Y/2 ) , L=lz, arrow=arrowaxes,p=black+.5bp); 

draw( shift(axesOrigin)*( O--.8(X+Z)/sqrt(8) ),  p=heavyred+1bp); 
draw( shift(axesOrigin)*( O--.8(Z-X)/sqrt(8) ) , p=heavygreen+1bp); 
draw( shift(axesOrigin)*( O--.8Y/2 ) , p=heavyblue+1bp); 