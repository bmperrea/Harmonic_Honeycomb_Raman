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
triple persp = (4,2.5,1.5);
real n = 4; //this is the ratio between how far away you thought you wanted the perspective and how far you actually want it. :)
currentprojection = perspective(n*persp,up=Z);

//Initial data

pen fontp = fontsize(14pt);
pen unitp = rgb(.5,0,.5)+1.5bp;
pen framep=gray(0.4)+.7bp;

material sphereCcolor = material(ambientpen= black, diffusepen=black, emissivepen=gray(.3));
//material(ambientpen= yellow, diffusepen=yellow, emissivepen=gray(.3));
material veccolor = material(diffusepen=red, ambientpen=red);
material framecolor = material(diffusepen=white,  ambientpen=white);
material sphereCcolor2 = material(ambientpen= yellow, diffusepen=yellow, emissivepen=gray(.2));
material sphereCcolor3 = material(ambientpen=gray(.7), diffusepen=gray(.7), emissivepen=gray(.4));
//rgb(255,182,193)

real cw = 2;
real diff = 0;

pen cxa = rgb(.8,.2,.2);
pen cxb = rgb(.2,.2,.8);
pen cya = rgb(.85,.5,0);
pen cyb = rgb(0,.5,.5);
pen cr = rgb(.4,.4,.4);
pen cb = rgb(0,0,0);


material cylcolorz = material(diffusepen=heavyblue+cw, ambientpen=deepblue);
material cylcolorxa = material(diffusepen=cxa+cw, ambientpen=deepred);//red
material cylcolorxb = material(diffusepen=cxb+cw, ambientpen=deepred);//blue
material cylcolorya = material(diffusepen=cya+cw, ambientpen=deepgreen);//orange
material cylcoloryb = material(diffusepen=cyb+cw, ambientpen=deepgreen);//green
material cylcolorg = material(diffusepen=heavygray+cw, ambientpen=darkgray);
material cylcolorr = material(diffusepen=cr+(cw+diff), ambientpen=deepblue);
material cylcolorb = material(diffusepen=cb+(cw-diff), ambientpen=black);

//pen cylcolorz = deepblue;
//pen cylcolory = deepgreen;
//pen cylcolorx = deepred;
//pen cylcolorg = darkgray;
//pen cylcolorr = deepblue;
//pen cylcolorb = gray;


real cylRadius = 0.1;
real sphereRadius = 0.05;
real sphereRadius2 = .14;
real frameCyl = .015;
real frameRadius = frameCyl;

//In these functions the smaller z value should come first

void drawRod(triple a, triple b, material cylcolor) {
	draw(a--b, cylcolor);
}

//void drawRod(triple a, triple b, material cylcolor) {
//  surface rod = extrude(scale(cylRadius)*unitcircle, axis=length(b-a)*Z);
//  triple orthovector = cross(Z, b-a);
//  if (length(orthovector) > .01) {
//    real angle = aCos(dot(Z, b-a) / length(b-a));
//    rod = rotate(angle, orthovector) * rod;
//  }
//  draw(shift(a)*rod, surfacepen=cylcolor);
//}

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
draw(shift(center)*scale3(sphereRadius2)*unitsphere, surfacepen=sphereCcolor2);
}

void drawCarbon3(triple center) {
draw(shift(center)*scale3(sphereRadius2)*unitsphere, surfacepen=sphereCcolor3);
}



//make a circle arc at p, with normal n and radius r
path3 circl(triple p, triple n, real r)
{
return shift(p)*rotate(-(180/pi)*acos(n.z/length(n)) ,cross(n,(0,0,1)))*scale3(r)*unitcircle3;
}



// lattice vectors 

triple dz = (0,0,1);
triple dxa = (-1,-sqrt(2),-1)/2;
triple dxb = (-1,sqrt(2),-1)/2;
triple dya = (1,sqrt(2),-1)/2;
triple dyb = (1,-sqrt(2),-1)/2;

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

//Unit Vectors
triple aa = Ca;
triple ab = Ba;
triple ac = Da;

aa = dxa - dya;
ab = dxb - dyb;
ac = 6*dz;



//Arcs for the gauge


real grad = .12;
pen gpen = gray(.15)+linewidth(1.0pt);

//triple point = bb-ab+dxb/2;
//draw(circl(point,dxb,grad),gpen);
//draw(circl(point+1.5*dz,dyb,grad),gpen);
triple point = O;




//Lattice sites

	triple bb = O;
	
	// rung with butt at bb
//	drawCarbon(bb);
//	triple bbxa = bb+dxa; drawCarbon(bbxa); drawRod(bbxa,bb,cylcolorx);
//	triple bbya = bb+dya; drawCarbon(bbya); drawRod(bbya,bb,cylcolory); 
//	triple bbr = bb+dz; drawCarbon(bbr); drawRod(bb,bbr,cylcolorr);  
//	triple bbrxa = bbr-dxa; drawCarbon(bbrxa); drawRod(bbr,bbrxa,cylcolorx);
//	triple bbrya = bbr-dya; drawCarbon(bbrya); drawRod(bbr,bbrya,cylcolory); 
	
//	bb = bb+ac;
	
	
		triple b=O;
		triple b2 = b + ab; triple b3 = b - ab;
		triple b02 = b + ac;
		triple b22 = b2 + ac; 
		triple b32 = b3 + ac;
	
	
	triple[] fluxps = {b2,b3,b22,b32};
	for(triple bb : fluxps) {
	
	point = bb-ab+dxa/2+ac/2+ab/2;
	draw(circl(point,dya,grad),gpen);
	draw(circl(point-4.5*dz,dxa,grad),gpen);
	
	draw(circl(point+aa,dya,grad),gpen);
	draw(circl(point+aa-4.5*dz,dxa,grad),gpen);
	}
	
		
		triple[] array = {b,b2,b3,b02,b22,b32};
		for(triple bb : array) {
	
	
	
	//H-1 B to B
	drawCarbon(bb);
	triple byb = bb+dyb; drawCarbon(byb); drawRod(byb,bb,cylcoloryb);
//	triple bxb = bb+dxb; drawCarbon(bxb); drawRod(bxb,bb,cylcolorxb);
	
	triple bybb = byb-dz; drawCarbon(bybb); drawRod(bybb,byb,cylcolorb);
//	triple bxbb = bxb-dz; drawCarbon(bxbb); drawRod(bxbb,bxb,cylcolorb);
	
	triple bybbya = bybb + dya; drawCarbon(bybbya); drawRod(bybb,bybbya,cylcolorya);
	triple bybbxa = bybb + dxa; drawCarbon(bybbxa); drawRod(bybb,bybbxa,cylcolorxa);
	 
	triple br = bb+dz; drawCarbon(br); drawRod(bb,br,cylcolorr);  
	triple brxb = br-dxb; drawCarbon(brxb); drawRod(br,brxb,cylcolorxb);
	triple bryb = brxb+dyb; drawCarbon(bryb); drawRod(bryb,brxb,cylcoloryb); 
	
	triple bxb = byb - dxb; drawCarbon(bxb); drawRod(bxb,byb,cylcolorxb);
	triple bxbr = bxb + dz; drawCarbon(bxbr); drawRod(bxb,bxbr,cylcolorr);
	
	triple brxbb = brxb + dz; drawCarbon(brxbb); drawRod(brxbb,brxb,cylcolorb);
	triple brxbbya = brxbb - dya; drawCarbon(brxbbya); drawRod(brxbbya,brxbb,cylcolorya);
	triple brxbbxa = brxbb - dxa; drawCarbon(brxbbxa); drawRod(brxbbxa,brxbb,cylcolorxa);
	
	
	triple bybbyar = bybbya - dz; drawCarbon(bybbyar); drawRod(bybbya,bybbyar,cylcolorr);
	triple bybbxar = bybbxa - dz; drawCarbon(bybbxar); drawRod(bybbxa,bybbxar,cylcolorr);
	
	
	
	}
			triple[] array2 = {b,b2,b3};
			for(triple bb : array2) {
	//b-type rungs
	triple btr = bb + 3*dz - dxb + dxa; drawCarbon(btr);
	triple btrz = btr + dz;  drawCarbon(btrz); drawRod(btr,btrz,cylcolorr);
	triple btr2 = btr-aa; drawCarbon(btr2);
	triple btr2z = btr2 + dz;  drawCarbon(btr2z); drawRod(btr2,btr2z,cylcolorr);
			}
	
////////////////////////////////////////////////////////////
//redefine all the colors to be dashed lines and then run the same stuff again moved over by aa.
pen dash = linetype(new real[] {2,2});
real w = .2;

material cylcolorz = material(diffusepen=heavyblue+cw+dash, ambientpen=deepblue);
material cylcolorxa = material(diffusepen=rgb(.8+w,.2+w,.2+w)+cw+dash, ambientpen=deepred);//red
material cylcolorxb = material(diffusepen=rgb(.2+w,.2+w,.8+w)+cw+dash, ambientpen=deepred);//blue
material cylcolorya = material(diffusepen=rgb(.85+w,.5+w,0+w)+cw+dash, ambientpen=deepgreen);//orange
material cylcoloryb = material(diffusepen=rgb(0+w,.5+w,.5+w)+cw+dash, ambientpen=deepgreen);//green
material cylcolorg = material(diffusepen=heavygray+cw+dash, ambientpen=darkgray);
material cylcolorr = material(diffusepen=rgb(.4+w,.4+w,.4+w)+(cw+diff)+dash, ambientpen=deepblue);
material cylcolorb = material(diffusepen=rgb(0+2*w,0+2*w,0+2*w)+(cw-diff)+dash, ambientpen=black);

b = b + aa;


		triple b2 = b + ab; triple b3 = b - ab;
		triple b02 = b + ac;
		triple b22 = b2 + ac; 
		triple b32 = b3 + ac;
		
		triple[] array = {b,b2,b3,b02,b22,b32};
		for(triple bb : array) {
		
		
		//H-1 B to B
		drawCarbon(bb);
		triple byb = bb+dyb; drawCarbon(byb); drawRod(byb,bb,cylcoloryb);
		//	triple bxb = bb+dxb; drawCarbon(bxb); drawRod(bxb,bb,cylcolorxb);
		
		triple bybb = byb-dz; drawCarbon(bybb); drawRod(bybb,byb,cylcolorb);
		//	triple bxbb = bxb-dz; drawCarbon(bxbb); drawRod(bxbb,bxb,cylcolorb);
		
		triple bybbya = bybb + dya; drawCarbon(bybbya); drawRod(bybb,bybbya,cylcolorya);
		triple bybbxa = bybb + dxa; drawCarbon(bybbxa); drawRod(bybb,bybbxa,cylcolorxa);
		
		triple br = bb+dz; drawCarbon(br); drawRod(bb,br,cylcolorr);  
		triple brxb = br-dxb; drawCarbon(brxb); drawRod(br,brxb,cylcolorxb);
		triple bryb = brxb+dyb; drawCarbon(bryb); drawRod(bryb,brxb,cylcoloryb); 
		
		triple bxb = byb - dxb; drawCarbon(bxb); drawRod(bxb,byb,cylcolorxb);
		triple bxbr = bxb + dz; drawCarbon(bxbr); drawRod(bxb,bxbr,cylcolorr);
		
		triple brxbb = brxb + dz; drawCarbon(brxbb); drawRod(brxbb,brxb,cylcolorb);
		triple brxbbya = brxbb - dya; drawCarbon(brxbbya); drawRod(brxbbya,brxbb,cylcolorya);
		triple brxbbxa = brxbb - dxa; drawCarbon(brxbbxa); drawRod(brxbbxa,brxbb,cylcolorxa);
		
		
		triple bybbyar = bybbya - dz; drawCarbon(bybbyar); drawRod(bybbya,bybbyar,cylcolorr);
		triple bybbxar = bybbxa - dz; drawCarbon(bybbxar); drawRod(bybbxa,bybbxar,cylcolorr);
		
		}
		triple[] array2 = {b,b2,b3};
		for(triple bb : array2) {
		//b-type rungs
		triple btr = bb + 3*dz - dxb + dxa; drawCarbon(btr);
		triple btrz = btr + dz;  drawCarbon(btrz); drawRod(btr,btrz,cylcolorr);
		triple btr2 = btr-aa; drawCarbon(btr2);
		triple btr2z = btr2 + dz;  drawCarbon(btr2z); drawRod(btr2,btr2z,cylcolorr);
		}




///////////////////////////////////////////////////////////
	
	
	
	// bridge with butt at bb
//	drawCarbon(bb);
//	triple bbxa = bb+dxa; drawCarbon(bbxa); drawRod(bbxa,bb,cylcolorxa);
//	triple bbya = bb+dya; drawCarbon(bbya); drawRod(bbya,bb,cylcolorya); 
//	triple bbr = bb+dz; drawCarbon(bbr); drawRod(bb,bbr,cylcolorb);  
//	triple bbrxb = bbr-dxb; drawCarbon(bbrxb); drawRod(bbr,bbrxb,cylcolorxb);
//	triple bbryb = bbr-dyb; drawCarbon(bbryb); drawRod(bbr,bbryb,cylcoloryb); 

//	triple bbxar = bbxa-dz; drawCarbon(bbxar); drawRod(bbxar,bbxa,cylcolorb);
//	triple bbyar = bbya-dz; drawCarbon(bbyar); drawRod(bbyar,bbya,cylcolorb);
//	}






//Alternative unit vectors that are useful
triple ad = -ab+aa;   //(goes along the ladder that a3 does not, doing so also with a positive b-component)
triple ae = ac - ad;  //(lies in a-direction)




//dimensions of the unit cell
real uX = 2*xpart(aa);
real uY = 2*ypart(ab);
real uZ = 2*zpart(ac);



//unit vectors

triple C=O-ac/2+aa/2+ab/2;

arrowbar3 arrowhead = Arrow3(size=150 sphereRadius);
margin3 themargin = Margin3(0,13 sphereRadius);

triple CC = C-ab;

draw(CC--(CC-aa),arrowhead,themargin,L=Label("$a_1$",0aa-1Y-1.6Z,p=fontp),p=unitp);
draw(CC--(CC-ab),arrowhead,themargin,L=Label("$a_2$",2ab-3.1Y+.4Z,p=fontp),p=unitp);
draw(CC--(CC+ac),arrowhead,themargin,L=Label("$a_3$",2ac-4.8Y-1.6Z,p=fontp),p=unitp);



//site names

drawCarbon2(CC);
drawCarbon3(CC+dz);
drawCarbon2(CC+dz-dxa);
drawCarbon3(CC+2*dz-dxa);
drawCarbon2(CC+2*dz-dxa-dyb);
drawCarbon3(CC+3*dz-dxa-dyb);
drawCarbon2(CC+3*dz-dxa-dyb-dxb);
drawCarbon3(CC+4*dz-dxa-dyb-dxb);





// Frame


draw(C--(C+(-uX,0,0)),framep);
draw(C--(C+(0,-uY,0)),framep);
draw(C--(C+(0,0,uZ)),framep);

draw((C+(-uX,-uY,0))--(C+(-uX,-uY,uZ)),framep);
draw((C+(-uX,-uY,uZ))--(C+(-uX,0,uZ)),framep);
draw((C+(-uX,-uY,uZ))--(C+(0,-uY,uZ)),framep);

draw((C+(-uX,-uY,0))--(C+(-uX,0,0)),framep);
draw((C+(-uX,-uY,0))--(C+(0,-uY,0)),framep);

draw((C+(-uX,0,0))--(C+(-uX,0,uZ)),framep);
draw((C+(-uX,0,uZ))--(C+(0,0,uZ)),framep);

draw((C+(0,-uY,0))--(C+(0,-uY,uZ)),framep);
draw((C+(0,-uY,uZ))--(C+(0,0,uZ)),framep);






//abc axes

real sc = 1.5;
real wid = 3;
real head = 0.5;

Label lx=Label("$a$",-1.2Z-.6X,p=fontp);
Label ly=Label("$b$",1.1Z+.3Y,p=fontp);
Label lz=Label("$c$",.95Y+.3Z,p=fontp);

triple axesOrigin = C;//-uX*X;
arrowbar3 arrowaxes = Arrow3(size=11);

draw( shift(axesOrigin)*( O--(-1)*X*sc*.8 ), L=lx ,arrow=arrowaxes,p=black+wid/1.5); 
draw( shift(axesOrigin)*( O--(-1)*Y*sc*.8) , L=ly, arrow=arrowaxes,p=black+wid/1.5); 
draw( shift(axesOrigin)*( O--Z*sc*.8 ) , L=lz, arrow=arrowaxes,p=black+wid/1.5); 


//xyz axes

Label lx=Label("$x$",.8Y+3.5Z-2X,p=fontp);
Label ly=Label("$y$",2.2Z+-0.6Y,p=fontp);
Label lz=Label("$z$",1.2Y+1.4Z,p=fontp);

axesOrigin = axesOrigin-uX*X-5Y+0X;
//arrowbar3 arrowaxes = Arrow3(size=9);

draw( shift(axesOrigin)*( O--(Z-X)*sc/sqrt(2) ), L=lx ,arrow=arrowaxes,p=black+head); 
draw( shift(axesOrigin)*( O--(X+Z)*sc/sqrt(2) ) , L=ly, arrow=arrowaxes,p=black+head); 
draw( shift(axesOrigin)*( O--.8Y*sc ) , L=lz, arrow=arrowaxes,p=black+head); 

triple of = -(2X-Y)/35; real sc2 = .65;
draw( shift(axesOrigin)*( (O-of)--(sc2*(X+Z)*sc/sqrt(2) -of) ),  p=cyb+wid); 
draw( shift(axesOrigin)*( (O-of)--(sc2*(Z-X)*sc/sqrt(2) -of) ) , p=cxb+wid); 
draw( shift(axesOrigin)*( (O+of)--(sc2*(X+Z)*sc/sqrt(2) +of) ),  p=cya+wid); 
draw( shift(axesOrigin)*( (O+of)--(sc2*(Z-X)*sc/sqrt(2) +of) ) , p=cxa+wid); 
draw( shift(axesOrigin)*( O--.5Y*sc ) , p=gray(.3)+1.2*wid); 