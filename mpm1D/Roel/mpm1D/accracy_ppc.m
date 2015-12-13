PPC = [4 6 8];
NEl = [512];
Err = zeros(3,length(NEl));

for i=1:3
   for j=1:length(NEl)
       Err(i,j) = vibrating_string_1D_mpm(NEl(j),PPC(i));
   end
end


semilogy(NEl, Err(1,:),'r','LineWidth',2)
hold on
semilogy(NEl, Err(2,:),'b','LineWidth',2)
hold on
semilogy(NEl, Err(3,:),'m','LineWidth',2)
xlabel('number of elements','Fontsize',12)
ylabel('RMS error [m]','Fontsize',12)
legend('4 PPC', '6 PPC', '8 PPC')

fileID = fopen('error4_large.txt','w');
fprintf(fileID,'%12s\r\n',Err(1,:));
fclose(fileID);
fileID = fopen('error6_large.txt','w');
fprintf(fileID,'%12s\r\n',Err(2,:));
fclose(fileID);
fileID = fopen('error8_large.txt','w');
fprintf(fileID,'%12s\r\n',Err(3,:));
fclose(fileID);