function F = fill_nan_field(F, X, Y)
% FILL_NAN_FIELD  Fill NaN points by linear interpolation from valid points,
% nearest-neighbour outside their convex hull (same idea as Section 5).
    bad = isnan(F);
    if ~any(bad(:)), return; end
    good = ~bad;
    F(bad) = griddata(X(good), Y(good), F(good), X(bad), Y(bad), 'linear');
    bad2 = isnan(F);
    if any(bad2(:))
        good2 = ~bad2;
        F(bad2) = griddata(X(good2), Y(good2), F(good2), X(bad2), Y(bad2), 'nearest');
    end
end
