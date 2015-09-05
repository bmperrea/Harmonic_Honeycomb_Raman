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
triple persp = (3*sqrt(2),3,-1);
real n = 4; //this is the ratio between how far away you thought you wanted the perspective and how far you actually want it. :)
currentprojection = perspective(n*persp,up=Z);

//Initial data

pen fontp = fontsize(24pt);
//pen unitp = rgb(.5,0,.5)+1.5bp;
//pen framep=gray(0.4)+.7bp;

//material sphereCcolor = material(ambientpen= black, diffusepen=black, emissivepen=gray(.3));
//material(ambientpen= yellow, diffusepen=yellow, emissivepen=gray(.3));
//material veccolor = material(diffusepen=red, ambientpen=red);
//material framecolor = material(diffusepen=white,  ambientpen=white);
//material sphereCcolor2 = material(ambientpen= yellow, diffusepen=yellow, emissivepen=gray(.2));
//material sphereCcolor3 = material(ambientpen=gray(.7), diffusepen=gray(.7), emissivepen=gray(.4));
//rgb(255,182,193)

real cw = 6;
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

//void drawFrame(triple a, triple b) {
//	surface rod = extrude(scale(frameCyl)*unitcircle, axis=length(b-a)*Z);
//	triple orthovector = cross(Z, b-a);
//	if (length(orthovector) > .01) {
//		real angle = aCos(dot(Z, b-a) / length(b-a));
//		rod = rotate(angle, orthovector) * rod;
//	}
//	draw(shift(a)*rod, surfacepen=framecolor);
//	draw(shift(b)*scale3(frameRadius)*unitsphere, surfacepen=framecolor);
//}

//void drawCarbon(triple center) {
//draw(shift(center)*scale3(sphereRadius)*unitsphere, surfacepen=sphereCcolor);
//}

//void drawCarbon2(triple center) {
//draw(shift(center)*scale3(sphereRadius2)*unitsphere, //surfacepen=sphereCcolor2);
//}

//void drawCarbon3(triple center) {
//draw(shift(center)*scale3(sphereRadius2)*unitsphere, surfacepen=sphereCcolor3);
//}



//make a circle arc at p, with normal n and radius r
//path3 circl(triple p, triple n, real r)
//{
//return shift(p)*rotate(-(180/pi)*acos(n.z/length(n)) //,cross(n,(0,0,1)))*scale3(r)*unitcircle3;
//}



// lattice vectors 

triple dz = (0,0,1);

triple dxb = (-1,sqrt(2),-1)/2;
triple dxa = dxb;//(-1,-sqrt(2),-1)/2;

triple dyb = (1,-sqrt(2),-1)/2;
triple dya = dyb;//(1,sqrt(2),-1)/2;


//Unit Vectors

triple aa = (-1,-sqrt(2),-1)/2 - (1,sqrt(2),-1)/2;
triple ab = dxb - dyb;
triple ac = 4*dz;

triple na = dz-dxa;
triple nb = dz-dya;


//labels
Label lx=Label("$x$",+.6nb,p=fontp);
Label ly=Label("$y$",-.6na,p=fontp);
Label lz=Label("$z$",+.5ab,p=fontp);

Label lx2=Label("$x$",-.6nb,p=fontp);
Label ly2=Label("$y$",+.6na,p=fontp);
Label lz2=Label("$z$",-.5ab,p=fontp);

triple bt = O; triple bt2 = bt - na;
triple btx = bt+dxa; draw(bt--btx, L=lx,p=cylcolorxb);
triple bt2y = bt2+dya; draw(bt2--bt2y, L=ly,p=cylcoloryb); 
triple bt2z = bt2+dz; draw(bt2--bt2z, L=lz,p=cylcolorr); 

triple bt3 = bt - nb;
triple bt3x = bt3+dxa; draw(bt3--bt3x, L=lx2,p=cylcolorxb);
triple bt3z = bt3+dz; draw(bt3--bt3z, L=lz2,p=cylcolorr);
triple bt3zy = bt3z-dya; draw(bt3z--bt3zy, L=ly2,p=cylcoloryb);

//labels
Label lz3=Label("$\sigma^z$",-6.7dz,p=fontp);
Label ly3=Label("$\sigma^y$",6.3dya,p=fontp);
Label lx3=Label("$\sigma^x$",-6.3dxa,p=fontp);
Label lz4=Label("$\sigma^z$",6.7dz,p=fontp);
Label ly4=Label("$\sigma^y$",-6.3dya,p=fontp);
Label lx4=Label("$\sigma^x$",6.3dxa,p=fontp);

triple bto = bt+dz; draw(bt--bto, L=lz3, p=cylcolorr);
triple btxo = btx-dya; draw(btx--btxo, L=ly3, p=cylcoloryb);
triple bt2o = bt2+dxa; draw(bt2--bt2o, L=lx3, p=cylcolorxb);
triple bt2yo = bt2y-dz; draw(bt2y--bt2yo, L=lz4, p=cylcolorr);
triple bt3o = bt3+dya; draw(bt3--bt3o, L=ly4, p=cylcoloryb);
triple bt3zo = bt3z-dxa; draw(bt3z--bt3zo, L=lx4, p=cylcolorxb);

