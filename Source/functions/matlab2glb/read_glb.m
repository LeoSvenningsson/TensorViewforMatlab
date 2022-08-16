function str = read_glb(fname)

    % focused on the meshed objects
    warning('read_glb is written for internal validation only and does not have full functionality')

    %% read

    fid = fopen(fname, 'r');
    str.binary = fread(fid);
    fclose(fid);
    
    %% read chunks
    
    to_uint16 = @(x) typecast(uint8(x), 'uint16');
    to_uint32 = @(x) typecast(uint8(x), 'uint32');
    to_single = @(x) typecast(uint8(x), 'single');    
    
    str.header.magic = char(str.binary(1:4)');
    str.header.version = to_uint32(str.binary(5:8));
    str.header.total_length = to_uint32(str.binary(9:12));
       
    c = 1;
    s = 13;
    
    while s < str.header.total_length

        chunks{c}.length = to_uint32(str.binary(s:s + 3));
        chunks{c}.type = to_uint32(str.binary(s + 4:s + 7));
        chunks{c}.data = str.binary(s + 8:s+chunks{c}.length + 7);
        s = s + chunks{c}.length + 8;
        c = c + 1;
        
    end
    
    %% parse chunks
    
    d = [];
    
    for c = 1:numel(chunks)        
        switch chunks{c}.type
            case 1313821514 % json
                str.json = jsondecode(char(chunks{c}.data'));
                to_cell = {'meshes' 'bufferViews' 'buffers' 'accessors'};
                % convert structures to cells so it's consistent every time
                for i = 1:numel(to_cell)
                    if isstruct(str.json.(to_cell{i}))
                        ne = numel(str.json.(to_cell{i}));
                        new_str = cell(ne, 1);
                        for j = 1:ne
                            new_str{j} = str.json.(to_cell{i})(j);
                        end
                        str.json.(to_cell{i}) = new_str;
                    end
                end
                d = [d c];
            case 5130562 % bin
                data = chunks{c}.data;
                d = [d c];
        end
    end
    
    chunks(d) = [];
    
    if ~isempty(chunks)
        str.other_chunks = chunks;
    end
    
    str = rmfield(str, 'binary');
    
    %% parse data    
    
    % external buffers    
    extbuff = cellfun(@(x)isfield(x, 'uri'), str.json.buffers);
    
    meshes = {};
           
    for i = 1:numel(str.json.meshes)
        
        % capture the name of the mesh
        if isfield(str.json.meshes{i}, 'name')
            meshes{i}.name = str.json.meshes{i}.name; 
        end
        
        % unwrap the mesh primitives
        if ~isfield(str.json.meshes{i}, 'primitives')
            str.json.meshes{i}.primitives{1} = str.json.meshes{i};
        elseif isstruct(str.json.meshes{i}.primitives)
            new_str = str.json.meshes{i}.primitives;
            str.json.meshes{i} = rmfield(str.json.meshes{i}, 'primitives');
            str.json.meshes{i}.primitives{1} = new_str;
        end 
        
        primitives = {};
        
        % go through each primitive
        for p = 1:numel(str.json.meshes{i}.primitives{1})

            p_raw = str.json.meshes{i}.primitives{p};

            % data types
            if isfield(p_raw, 'attributes')
                p_data = p_raw.attributes;
            else                
                p_data = struct;
            end
            if isfield(p_raw, 'indices')
                p_data.indices = p_raw.indices;
            end
            
            fn = fieldnames(p_data);
            if isempty(fn)
                continue
            end
            
            % check for external data - not parsing for now 
            if any(extbuff)            
                accessors = cellfun(@(x)(p_data.(x)), fn) + 1;
                bufferViews = arrayfun(@(x) str.json.accessors{x}.bufferView, accessors) + 1;
                buffers = arrayfun(@(x) str.json.bufferViews{x}.buffer, bufferViews) + 1;
                if any(extbuff(buffers))
                    warning('Primitive uses external data - not parsed'); 
                    continue
                end
            end

            % parse each type    
            for j = 1:numel(fn)

                % choose the correct accessor and note down data name
                tag = fn{j};
                a = p_data.(fn{j}) + 1;                
                
                switch tag
                    case 'indices'
                        tag = 'F';
                    case 'POSITION'
                        tag = 'V';
                end

                % get the pointers
                bw  = str.json.accessors{a}.bufferView + 1;          % bufferView
                bwo = str.json.bufferViews{bw}.byteOffset;           % bufferView offset
                bwl = str.json.bufferViews{bw}.byteLength;           % bufferView length               
                acc = str.json.accessors{a}.count;                   % accessors count 
                if isfield(str.json.accessors{a}, 'byteOffset')      % accessors offset
                    aco = str.json.accessors{a}.byteOffset;              
                else
                    aco = 0;
                end
                
                % work out number of components per element
                switch str.json.accessors{a}.type                    % accessors type
                    case 'SCALAR'
                        cpe = 1;                        
                    case 'VEC2'
                        cpe = 2;
                    case 'VEC3'
                        cpe = 3;
                    case 'VEC4'
                        cpe = 4;                       
                    otherwise
                        error('accessor type not handled at present')
                end
                
                % work out number of bytes per component
                ctp = str.json.accessors{a}.componentType;           % accessors component type
                nbt = 2 ^ (floor((mod(ctp, 5120) - 1) / 2));  
                
                % work out the stride
                if isfield(str.json.bufferViews{bw}, 'byteStride')                  
                    bws = str.json.bufferViews{bw}.byteStride;       % bufferView stride
                else
                    bws = cpe * nbt;
                end
                
                % extract relevant bytes
                bv_data = reshape(data(bwo + 1:bwo + bwl), bws, []);
                [r, c] = ind2sub([bws bwl], aco + 1);
                ac_data = reshape(bv_data(r:r + cpe * nbt - 1, c:c + acc - 1), [], 1);

                % parse as appropriate
                switch str.json.accessors{a}.componentType
                    case 5121
                        parsed_data = ac_data;
                    case 5123
                        parsed_data = to_uint16(ac_data);
                    case 5125
                        parsed_data = to_uint32(ac_data);
                    case 5126
                        parsed_data = to_single(ac_data);
                    otherwise
                        error('add other options')
                end            

                if strcmp(tag, 'F')
                    parsed_data = parsed_data + 1; % add 1 to indices - matlab compatibility
                    parsed_data = reshape(parsed_data', 3, [])';                
                else
                    parsed_data = reshape(parsed_data', cpe, [])';
                end
                
                parsed_data = double(parsed_data); % for matlab

                primitives{p}.(tag) = parsed_data;

            end          
        end 
        
        % unwrap if only one primitive is present
        if numel(primitives) > 1
            meshes{i}.primitives = primitives;
        else
            fn = fieldnames(primitives{1});
            for f = 1:numel(fn)
                meshes{i}.(fn{f}) = primitives{1}.(fn{f});
            end
        end
    end
    
    if numel(meshes) == 1
        meshes = meshes{1};
    end
    
    str.mesh = meshes;    
    
end