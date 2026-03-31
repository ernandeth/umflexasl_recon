kv=load('kviews27648.txt');

unit = [1 0 1]
for n=1:size(kv,1)
    tmp = kv(n, end-8:end);
    R = reshape(tmp',3,3)';
    s = unit*R';
    plot3(s(1),s(2), s(3), 'o');
    hold on
    drawnow
end
hold off
