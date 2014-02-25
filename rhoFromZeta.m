function [rho3d] = rhoFromZeta( zeta3d, z, rho_bar)
	rho3d = zeros(size(zeta3d));
	for m=1:size(rho3d,1)
		for n=1:size(rho3d,2)
			coordinate = squeeze(zeta3d(m,n,:))+z;
			negIndex = find(diff(coordinate)<=0);
			if (length(negIndex) ~= 0)
				for j=(min(negIndex)+1):(max(negIndex)+1)
					if (coordinate(j)<=coordinate(j-1))
						coordinate(j)=coordinate(j-1)+0.001;
					end
				end
			end
			rho3d(m,n,:) = interp1( coordinate, rho_bar, z );
		end
	end
	rho3d(find(isnan(rho3d)))=min(rho_bar);
end