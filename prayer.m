% This is a script I wrote for calculating how close I was to 99 prayer on RuneScape.
% it took the XP left to 99 prayer and converted it into a % from my start point at ~200-300k til 97 to 99
% along with a number of trips to the bank back to my gilded altar left and the number of dragon bones I
% need to offer at my gilded altar. 
function sol = prayer(t)
  [(1 - t / (252*1.06*8400))*100, t/(252*1.06*55), t/(252*1.06)]
