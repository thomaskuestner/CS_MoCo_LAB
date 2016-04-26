function Plot3DData(S,gridSize,optionString)

for i=1:gridSize
  for j=1:gridSize
    for k=1:gridSize
      if S(i,j,k) >=.5
        plot3(i,j,k,optionString);
      end
    end
  end
end
      