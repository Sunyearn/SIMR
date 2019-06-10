clear;
%%
[model,block,per] = creat_model_commemi_3d_1();
point = struct; % initiate the struct "point", which contains position, value of physical field, etc.
line_x = struct; % initiate the struct "line_x", which contains position, value of physical field, etc.(line in x direction)
line_y = struct; % initiate the struct "line_y", which contains position, value of physical field, etc.(line in y direction)
for ilevel = 1:(model(1).level(end))
    for iblock = 1:(block(ilevel).flag(1,2))
        [point,line_x,line_y] = creat_measure_point(block,ilevel,iblock,point,line_x,line_y);
    end
end
Nper = length(per); % the number of vector of period

tic
skip = [0,0,0];
for iper = 1:Nper
    for ilevel = 1:(model(1).level(2))
        if ilevel == 1
            [field] = mtfwd3d_initial(block,ilevel,1,per(iper),'bicg',1000,1e-10); %'bicg': biconjugate gradient stabilized method
            [model,block] = fill_in_field2model(model,block,ilevel,1,iper,field,skip);
            continue;
        else
            for iblock = 1 : (block(ilevel).flag(1,2))
                [model,block] = extend_e(model,block,ilevel,iblock,iper);
                [field]= mtfwd3d3(block,ilevel,iblock,iper,per(iper),'bicg',1000,1e-10);
                [model,block] = fill_in_field2model(model,block,ilevel,iblock,iper,field,skip);
            end
        end
    end
end
toc
save('model.mat','model');
save('block.mat','block');
% function plot_model_or_block(field3D,ilevel,iblock,iper,per,ipol,ilayer,type,mode,view,padding,model_or_block)
% plot_model_or_block(block,1,1,1,per,1,1,'E','X','Z',padding,'BLOCK')
% plot_model_or_block(block,2,1,1,per,1,1,'E','X','Z',padding,'BLOCK')
% plot_model_or_block(block,2,2,1,per,1,1,'E','X','Z',padding,'BLOCK')

for ilevel = 1:(model(1).level(2))
    for iblock = 1:(block(ilevel).flag(1,2))
        [point] = cal_impedance(block,point,ilevel,iblock,per);
        [point] = cal_rho(point,ilevel,iblock,per);
        [line_x] = cal_impedance(block,line_x,ilevel,iblock,per);
        [line_x] = cal_rho(line_x,ilevel,iblock,per);
        [line_y] = cal_impedance(block,line_y,ilevel,iblock,per);
        [line_y] = cal_rho(line_y,ilevel,iblock,per);
    end
end


for iper = 1:length(per)
    for ilevel = 1:(model(1).level(2))
        for iblock = 1:(block(ilevel).flag(1,2))
            plot_rho_theta(line_x,ilevel,iblock,iper,per,'LINE','X');
            plot_rho_theta(line_y,ilevel,iblock,iper,per,'LINE','Y');
        end
    end
end