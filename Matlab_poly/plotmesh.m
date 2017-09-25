function [] = plotmesh(mesh)
%[] = plotmesh(mesh)
%   Plots the mesh.

n_step = ceil(mesh.K^(1/3));

figure
for ie = 1:mesh.Ntet
    conv = dec2base(mesh.E2P(ie)-1, n_step,3);
    s1 = str2double(conv(1));
    s2 = str2double(conv(3));
    s3 = str2double(conv(2));
    c = [s1 s2 s3]/n_step;
    for s = 1:4
        for t = s+1:4
            plot3(mesh.VX(mesh.E2V(ie, [s t])), mesh.VY(mesh.E2V(ie, [s t])), mesh.VZ(mesh.E2V(ie, [s t])), 'Color', c);
            hold on
        end
    end
end
hold off

grid on

end

