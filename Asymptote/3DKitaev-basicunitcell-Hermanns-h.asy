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

pen colorz = blue;
pen colory = green;
pen colorx = red;
material cylcolorz = material(diffusepen=blue, ambientpen=deepblue);
material cylcolory = material(diffusepen=green, ambientpen=deepgreen);
material cylcolorx = material(diffusepen=red, ambientpen=deepred);
material cylcolorg = material(diffusepen=lightgray, ambientpen=darkgray);

real cylRadius = 0.033;
real sphereRadius = 0.2;
real frameCyl = .015;
real frameRadius = frameCyl;

//In these functions the smaller z value should come first

arrowbar3 arrowJ = Arrow3(size=15);

void drawRod(triple a, triple b, material cylcolor) {
  surface rod = extrude(scale(cylRadius)*unitcircle, axis=length(b-a)*Z);
  triple orthovector = cross(Z, b-a);
  if (length(orthovector) > .01) {
    real angle = aCos(dot(Z, b-a) / length(b-a));
    rod = rotate(angle, orthovector) * rod;
  }
  draw(shift(a)*rod, surfacepen=cylcolor);
  draw(a--(b+0.05*(a-b)), arrow=arrowJ, p=white+1bp);
}

pen dotted=linetype(new real[] {8,8});
arrowbar3 arrowk = Arrow3(size=15);

void drawRod2(triple a, triple b, pen col) {
//surface rod = extrude(scale(cylRadius)*unitcircle, axis=length(b-a)*Z);
//triple orthovector = cross(Z, b-a);
////if (length(orthovector) > .01) {
real angle = aCos(dot(Z, b-a) / length(b-a));
//rod = rotate(angle, orthovector) * rod;
//}
//draw(shift(a)*rod, surfacepen=cylcolor);
draw(a--(b+0.05*(a-b)), arrow=arrowk, p=col+1.5bp+dotted);
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

// lattice vectors

triple dz = (0,0,sqrt(2));
triple dxa = (-1/sqrt(3),-sqrt(2/3),-1);
triple dxb = (-1/sqrt(3),sqrt(2/3),-1);
triple dya = (1/sqrt(3),sqrt(2/3),-1);
triple dyb = (1/sqrt(3),-sqrt(2/3),-1);

//Lattice sites

triple Aa = O;
triple Ab = dz;
triple Ac = Ab - dxa;
triple Ad = Ac + dz;
triple Ba = Ad - dyb;
triple Ca = Ad - dxb; 
triple Db = Ac + dya;
triple Da = Db - dz;
triple Ea = (2 xpart(Da),0,0);
triple Eb = Ea + dz;
triple Fa = Ba + Ea;

label("$1$",Aa,persp-Aa/n,p=fontp);
label("$2$",Ab, persp-Ab/n,p=fontp);
label("$3$",Ac, persp-Ac/n,p=fontp);
label("$4$",Ad, persp-Ad/n,p=fontp);

//Unit Vectors
triple aa = -dxb+dyb;
triple ab = -dxa+dya;
triple ac = -dxa-dxb+2dz;

arrowbar3 arrowhead = Arrow3(size=80 sphereRadius);
margin3 themargin = Margin3(0,13 sphereRadius);

draw(Da--(Da+aa),arrowhead,themargin,L=Label("$a_1$",2aa+4Y-2Z,p=fontp),p=unitp);
draw(Ab--(Ab+ab),arrowhead,themargin,L=Label("$a_2$",2ab-3Y+2.2Z,p=fontp),p=unitp);
draw(Aa--ac,arrowhead,themargin,L=Label("$a_3$",2ac-4Y-.6Z,p=fontp),p=unitp);

//dimensions of the unit cell
real uX = xpart(Ea);
real uY = ypart(Ba);
real uZ = zpart(Ba);


//Bonds

//blue (z)
drawRod(Aa,Ab,cylcolorz);
drawRod(Ac,Ad,cylcolorz);

//repeats (gray)
drawRod(Da,Db,cylcolorg);
drawRod(Ea,Eb,cylcolorg);

//red (x)
drawRod(Ac,Ab,cylcolorx);
drawRod(Ca,Ad,cylcolorx);

//green (y)
drawRod(Ac,Db,cylcolory);
drawRod(Ba,Ad,cylcolory);


	// The second neighbor arrows
drawRod2(Aa,Ac,colory);
drawRod2(Ab,Ad,colorz);
drawRod2(Ac,Da,colorx);
drawRod2(Ac,Ca,colory);
drawRod2(Ba,Ac,colorx);
drawRod2(Ca,Ba,colorz);
//drawRod2(Ea,Da,colorz);
drawRod2(Db,Ad,colorx);
drawRod2(Ab,Db,colorz);




//The labelled sites (yellow)
drawCarbon(Aa);
drawCarbon(Ab);
drawCarbon(Ac);
drawCarbon(Ad);

// The other ones (white)
drawCarbon2(Ba);
drawCarbon2(Ca);
drawCarbon2(Da);
drawCarbon2(Db);
drawCarbon2(Ea);
drawCarbon2(Eb);
drawCarbon2(Fa);



// Frame

draw((0,0,0)--(uX,0,0),framep);
draw((0,0,0)--(0,uY,0),framep);
draw((0,0,0)--(0,0,uZ),framep);

draw((uX,uY,0)--(uX,uY,uZ),framep);
draw((uX,uY,uZ)--(uX,0,uZ),framep);
draw((uX,uY,uZ)--(0,uY,uZ),framep);

draw((uX,uY,0)--(uX,0,0),framep);
draw((uX,uY,0)--(0,uY,0),framep);

draw((uX,0,0)--(uX,0,uZ),framep);
draw((uX,0,uZ)--(0,0,uZ),framep);

draw((0,uY,0)--(0,uY,uZ),framep);
draw((0,uY,uZ)--(0,0,uZ),framep);


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

Label lx=Label("$x$",.5Y+2.8Z,p=fontp);
Label ly=Label("$y$",1.7Z+.6Y,p=fontp);
Label lz=Label("$z$",1.2Y+Z,p=fontp);

triple axesOrigin = (-1.5,2,2.5);
arrowbar3 arrowaxes = Arrow3(size=7);

draw( shift(axesOrigin)*( O--(Z-X)/sqrt(8) ), L=lx ,arrow=arrowaxes,p=black+.5bp); 
draw( shift(axesOrigin)*( O--(X+Z)/sqrt(8) ) , L=ly, arrow=arrowaxes,p=black+.5bp); 
draw( shift(axesOrigin)*( O--Y/2 ) , L=lz, arrow=arrowaxes,p=black+.5bp); 

draw( shift(axesOrigin)*( O--.8(X+Z)/sqrt(8) ),  p=heavyred+1bp); 
draw( shift(axesOrigin)*( O--.8(Z-X)/sqrt(8) ) , p=heavygreen+1bp); 
draw( shift(axesOrigin)*( O--.8Y/2 ) , p=heavyblue+1bp); 