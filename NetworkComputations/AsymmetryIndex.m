function Index = AsymmetryIndex(Diff)
    Flipped_Diffs = Diff-flip(Diff);
    Index = norm(Flipped_Diffs);
    if sum(Flipped_Diffs(1:length(Diff)/2)) > sum(Flipped_Diffs((length(Diff)/2)+1:end))
        Index = -Index;
    elseif sum(Flipped_Diffs(1:length(Diff)/2)) < sum(Flipped_Diffs((length(Diff)/2)+1:end))
        Index = Index;
    else
        Index = Index;
        display("No directionality to asymmetry")
    end

end

