function [] = matrix_to_txt(mat,file_name,header)

if nargin < 1 
	file_name = 'mat_to_txt.txt';
end

T = size(mat,1);
N = size(mat,2);

% open a file for writing
fid = fopen(file_name, 'w');

if nargin > 2 % header found 
	for j = 1:size(header,1)
		fprintf(fid,	'%s\t',  strtrim(header(j,:)) );
	end;
	fprintf(fid, '\n'); % on ferme la ligne
end

% now writing the dataset
for t = 1:T
	for j = 1:N
		fprintf(fid,	'%f\t', mat(t,j) );
	end;
	fprintf(fid, '\n'); % on ferme la ligne
end;


fclose(fid);



