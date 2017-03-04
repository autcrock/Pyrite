%% ****************************************************************%
%% This program simulates the formation of pyrite crystals on a    %
%% 3 dimensional cubic lattice. The uncombined components diffuse  %
%% on this lattice. Simple random diffusion is provided for,       %
%% with reaction being immediate if the geometry is correct.       %
%% ****************************************************************%
%%
%% declare some things                              mt 20/5/85
%%
 
3 constant print-width
 
variable number-of-ions       %%
variable water+iron           %%
variable water/iron           %% various amounts for initialization
variable water/sulphide       %%
 
0 constant water              %%
1 constant iron               %% various ions in the solution
2 constant sulphide           %%
 
fvalue xyz                 %% current position of lattice pointer
 
value suitable              %% contents of new cell
value okvar                 %% the value of the reaction candidate + 10
value number-of-looks       %% how many times should random-look look
 
defer init        %% initialisation of lattice
defer set-cell    %% used by 3d!
defer action      %% used by 3dact
 
%%
%%
%% random number words
%%
%%
 
value seed here to seed
 
: random  seed 31421 *  6927 + dup to seed ;
 
: choose     (s u1  --- 0<=u2<u1 )
      random  um* swap drop ;
 
%%
%% lattice declaration                                 mt        
%%
 
20 dup dup
constant xmax
constant ymax
constant zmax
 
>memory
dim xmax ymax zmax byte array lattice
>dictionary
 
%%
%% Concise array access words
%%
 
: cell@     (s x y z --- contents )
      lattice pos c@ ;
 
: A$Kz%% 3d-array display words                              mt        
%%
 
: 1d.    (s x y  --- )
     cr zmax 0 do
         2dup i cell@  print-width u.r loop 2drop ;
 
: 2d.     (s x --- )
     cr ymax 0 do
         dup i 1d. loop drop ;
 
: 3d.     (s --- )
     cr xmax 0 do
         cr cr ." Slice number " i 1+ u.
         i 2d. loop ;
 
: slice   (  n ---) %% displays a single slice of the lattice
     cr cr ." Slice number " dup u. 1- cr 2d. ;
 
%%
%% 3d-array initialisation words                                 29/4/85   
%%
 
: 1d!    (s x y  --- )
      zmax 0 do
        2dup i set-cell
      loop 2drop ;
 
: 2d!     (s x --- )
      ymax 0 do
         dup i 1d! loop drop ;
 
: 3d!     (s --- )
      xmax 0 do
         i 2d! loop ;
 
%%
%% 3d-array processing words                                 mt 20/6/85
%% action can be set to anything
%%
 
: 1dact    (s x y  --- )
      zmax 2- 2 do
        2dup i action
      loop 2drop ;
 
: 2dact     (s x --- )
      ymax 2- 2 do
         dup i 1dact loop drop ;
 
: 3dact     (s --- )
      xmax 2- 2 do
         i 2dact loop ;
 
%%
%% random lattice filler                               mt 20/5/85
%%
%% Fills the lattice with (spacially) randomly distributed ions
%% The ratios of ions are held at the input values
%%
 
: coord-roll   3 roll 3 roll 3 roll ;
 
: random-ion-put  (s x y z --- )
    number-of-ions @ choose dup
          100 <
              if drop iron else
          water+iron @ <
              if water else sulphide then
              then
    coord-roll cell! ;
 
%%
%% words to increment/decrement xyz coords; mt 29/4/85
%%
 
: xinc   rot 1+ -rot ;
: xdec   rot 1- -rot  ;
: yinc   swap 1+ swap ;
: ydec   swap 1- swap ;
: zinc   1+ ;
: zdec   1- ;
 
%%
%% words to get contents surrounding the lattice pointer
%%
 
: xget  (s x y z --- next )
   xinc  cell@ ;
: yget  (s x y z --- next )
    yinc   cell@ ;
: zget  (s x y z --- next )
    zinc  cell@ ;
 
: -xget  (s x y z --- last )
   xdec cell@ ;
: -yget  (s x y z --- last )
   ydec  cell@ ;
: -zget  (s x y z --- last )
    zdec  cell@ ;
 
: xyget  (s x y z --- next )
     xinc yinc  cell@ ;
: -xyget      xdec yinc cell@ ;
: x-yget     xinc ydec cell@  ;
 
: xzget  (s x y z --- next )
       xinc zinc  cell@ ;
: -xzget   xdec zinc cell@ ;
: x-zget     xinc zdec cell@ ;
 
: yzget  (s x y z --- next )
    yinc zinc  cell@ ;
: -yzget     ydec zinc cell@ ;
: y-zget     yinc zdec cell@ ;
 
: -x-yget  (s x y z --- last )
   xdec ydec  cell@ ;
: -x-zget  (s x y z --- last )
    xdec zdec  cell@ ;
: -y-zget  (s x y z --- last )
    ydec zdec  cell@ ;
 
%%
%% prelim words for rule actions                       mt 29/4/85
%%
%% First, access the value at the lattice pointer
%%
 
: extract  xyz cell@ ;   (s  --- n )
 
%%
%% Now decide in which random direction to diffuse,
%% relative to the pointer.
%%
 
: +or-   2 choose negate ?negate ;   (s n --- +/-n )
 
: xchange    (s x y z  --- newx y x )
    rot  2 choose +or- + -rot ;
: ychange   (s x y z --- x newy z )
    swap 2 choose +or- + swap ;
: zchange     (s x y z --- x y newz )
         2 choose +or- + ;
 
%%
%% Check neighbours for a geometrically suitable bonding position,
%% and act accordingly. ie Either diffuse or react.
%%                             mt 29/4/85
%%
 
: cell-swap      (s x y z x1 y1 z1 --- )
     fdup cell@ >r   fswap fdup cell@
     r> swap >r coord-roll cell! r> coord-roll cell! ;
 
: destination     xyz xchange     ychange     zchange ;
 
: diffuse   (s  --- )  %% exchange random neighbour
    destination fdup cell@ 9 > if
         fdrop destination fdup cell@ 9 > if
              fdrop else xyz cell-swap then  else
    xyz cell-swap then ;
 
%%
%% Ions held in a crystal are indicated by a value of 10+original value
%%
 
: ?10+     (s   n1 --- n1 or n1+10 )
    dup 9 < if 10 + then ;
 
%%
%% react in various directions                         mt        
%%
 
: react     (s   ---  )
    xyz  cell@ ?10+   xyz cell! ;
: xreact      (s   --- )
    xyz xinc cell@ ?10+ xyz xinc cell! ;
: yreact      (s   --- )
    xyz yinc cell@ ?10+ xyz yinc cell! ;
: zreact      (s   --- )
    xyz zinc cell@ ?10+ xyz zinc cell! ;
: -xreact      (s   --- )
    xyz xdec cell@ ?10+ xyz xdec cell! ;
: -yreact      (s   --- )
    xyz ydec cell@ ?10+ xyz ydec cell! ;
: -zreact      (s   --- )
    xyz zdec cell@ ?10+ xyz zdec cell! ;
 
%%
%% check if surrounding ions prohibit a bond           mt 30/4/85
%% leave 0 if can't react, -1 otherwise
%%
 
: ?okr       (s x y z --- 0 or 1 )
    fdup cell@ 10 + to okvar
    fdup  xget okvar = if fdrop 0 else
    fdup  yget okvar = if fdrop 0 else
    fdup  zget okvar = if fdrop 0 else
    fdup -xget okvar = if fdrop 0 else
    fdup -yget okvar = if fdrop 0 else
         -zget okvar = if       0 else
    -1 then then then then then then ;
 
%%
%% check neighbours, act                               mt 30/4/85
%%
 
: look-around    (s    ---   )
   suitable xyz xget  = if xyz fdup xinc ?okr coord-roll
        ?okr and if react  xreact then else
   suitable xyz yget  = if xyz fdup yinc ?okr coord-roll
        ?okr and if react  yreact then else
   suitable xyz zget  = if xyz fdup zinc ?okr coord-roll
        ?okr and if react  zreact then else
   suitable xyz -xget = if xyz fdup xdec ?okr coord-roll
        ?okr and if react -xreact then else
   suitable xyz -yget = if xyz fdup ydec ?okr coord-roll
        ?okr and if react -yreact then else
   suitable xyz -zget = if xyz fdup zdec ?okr coord-roll
        ?okr and if react -zreact then else
   diffuse then then then then then then ;
 
%%
%% carry out the action for one lattice point           mt        
%%
 
:  random-diffuse-immediate-react       (s   x y z --- )
  to xyz     extract
  dup 9    > if drop                                   else
  dup iron = if drop sulphide to suitable look-around    else
  sulphide = if      iron     to suitable look-around    else
  diffuse
  then then then ;
 
%%
%% process the lattice                                 mt 21/5/85
%%
 
: simple (s   --- )
    ['] random-diffuse-immediate-react is action ;
 
: do-simple   (s --- )
    simple
    3dact ;
 
%% ***************************************************************%
 
 
%%
%% Random check-around, to avoid any artificial patterns
%%
%% make-delta-xyz sets up the changes in the three coords
%% add-delta-xyz adds the changes to the coord.
%% random-get gets any-old-ion
%% random-look looks around till it runs out of time or reacts
%%
 
: make-delta-xyz      ( --- n1 n2 n3 , all 0 or +/-1)
    2 choose +or-
    2 choose +or-
    2 choose +or- ;
 
: add-delta-xyz       ( x y z n1 n2 n3 --- x' y' z')
    >r >r >r          ( x y z )
    rot r> +          ( y z x')
    rot r> +          ( z x' y')
    rot r> +  ;       ( x' y' z')
 
: random-get          ( --- x' y' z' ion )
    xyz make-delta-xyz add-delta-xyz fdup c@ ;
 
: random-look          ( --- )
  number-of-looks 1 max 0 do
    random-get suitable = if
        fdup ?okr xyz ?okr and if
            fdup c@ ?10+ c!
            xyz fdup c@ ?10+ c!
            then then loop ;
 
%%
%% carry out the action for one lattice point           mt        
%%
 
:  rand-dif-imm-react-rand-look       (s   x y z --- )
  to xyz     extract
  dup 9    > if drop                                   else
  dup iron = if drop sulphide to suitable random-look    else
  sulphide = if      iron     to suitable random-look    else
  diffuse
  then then then ;
 
 
%%
%% process the lattice                                 mt 21/5/85
%%
 
: simple1 (s   --- )
    ['] rand-dif-imm-react-rand-look is action ;
 
: do-simple1   (s --- )
    simple1
    3dact ;
 
%% **************************************************************%
 
%%
%% get the component ratios from the console      mt 22/5/85
%%
 
: getratio   (s --- )
    cr ." Input the molar ratio of water:iron. "
    sinp dup water/iron !
    100 * dup
    100 + dup water+iron !  swap
    cr ." Input the molar ratio of water:sulphide. "
    sinp dup water/sulphide !
    / + number-of-ions !  ;
 
%%
%% Initialise the lattice initialisation word
%%
 
: ions   ['] random-ion-put is set-cell ;
 
%%
%% Set up the lattice using the random model
%%
 
: set-random
    getratio ions  3d! ;
 
%%
%% multiple cycles                   mt 22/5/85
%%
 
: cycles    (s n --- )
      set-random  3d.
      1 max 0 do do-simple 3d. loop ;
 
: cycle cycles ;
