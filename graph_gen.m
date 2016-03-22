%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Graph Generator v1.0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%											%
% 				    University of York					%
% 				     Graceful Project					%
% 				  Created by Nizar Dahir				%
% 					York, 2016					%
%											%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	clear
	% Set the seed to repeat results
	seed = 15;
	rand ("seed", seed)
	randn("seed", seed)

	nodes =100;
	
	min_actual_ranks = 2*sqrt(nodes);
	minrank = floor(min_actual_ranks*1.2) ;
	maxrank = floor(min_actual_ranks*1.5);

	maxload = 10;	

	% Number of graphs		
	N = 20; 
	
	% Probability of connection
	prob=0.1;

	for g=1:N
		valid = 0;
		while (~valid)
			conn = zeros(nodes,nodes);

			ranks = randi([minrank maxrank],1);
			ranks = floor(ranks + 0.5);

			% mean of node binning
			mu_rank =  ranks/2;

			% variance of node binning
			sig_rank =  ranks/3;

			while (1)
				for i=1:nodes
					nranks(i) =  limitednorm(mu_rank,sig_rank, 0.5, ranks + 0.5);
%					nranks(i) =  randi([1 ranks], 1);
				end
	
				nranks=  floor(nranks + 0.5);

				nranks=  sort(nranks);	

				if (length(unique(nranks)) >= min_actual_ranks) break; end
		
			end

			for s=1 : nodes -1
				for d= s + 1 :nodes
					conn(s,d) =  nranks(d) - nranks(s);
				end	
			end


			for s=1 : nodes - 1
				for d= s + 1 :nodes
					if (conn(s,d) > 0) conn(s,d) = prob^conn(s,d); end
				end
			end

			conn = triu(conn,1);
		
			conn = (rand(nodes) < conn);

			% nonempty ranks
			aranks=unique(nranks);

			ideg = sum(conn); 
			odeg = sum(conn');

			% repair unconnected nodes
			for i=1:nodes
				rnode = nranks(i);
				Fi = find(aranks == rnode);
				if (~ideg(i) && (nranks(i) > aranks(1)))			
					Fi = find(aranks == rnode);
					prev_rank = aranks(Fi - 1);
					prev_nodes = find(nranks== prev_rank);
					[Av, Ai] = min(odeg(prev_nodes));	
					conn(prev_nodes(Ai),i) = 1;
				end
				if (~odeg(i) && (rnode < aranks(end)))			
					Fi = find(aranks == rnode);
					next_rank = aranks(Fi + 1);
					next_nodes = find(nranks == next_rank);
					[Av, Ai] = min(ideg(next_nodes));	
					conn(i, next_nodes(Ai)) = 1;
				end

			end
	



			% check if the graph has any unconnected subgraphs 
			n=[];
			V=[];
			V0=1;		
			while (1)
				for i=1:length(V0)
				 	n = [n find(conn(V0(i),:)) find(conn(:,V0(i)))'];
				end

				Vn = setdiff(n, V);
				V =[V Vn];
				if (isempty(Vn)) break; end
				V0=Vn;	
			end

			if (length(V) == nodes)
				valid = 1;
			else
				disp("invalid !");
			end

		end  % while (~valid)

			 
		g;
		conn ;

		
		tdeg = sum(conn) + sum(conn');

		
		% Write gv file 
		fname = sprintf("graphs/dag%04i",g);
		fid=fopen( [fname '.gv'] , "w");

		str = ["digraph {\n"];
		str = [str "  splines=true;\n\r"];
		str=[str "node [margin=0 fontname=arial fontcolor=black fontsize=12 shape=circle "]; 
		str=[str "width=0.5 fixedsize=true style=filled fillcolor=powderblue]\n\r"];

		% Nodes
		for i=1:nodes
			str = [str sprintf("  %i [label=\"P%i\"]\n",i,i)];
		end
	
		str = [str sprintf("rankdir=LR\n\r")];
		str = [str sprintf("edge [margin=0 fontname=arial fontcolor=black fontsize=12]\n\r")];

		[s d] = find(conn==1);

		edges = length(s); 

		% Edges 
		loads = randi([1 maxload], edges, 1);

		for i=1:length(s) 
			str = [str sprintf("	%i -> %i [label=\"%i\"]\n",s(i),d(i), loads(i))];

		end
		
		% Ranks
		for i=aranks
			same = find(nranks==i);

			str = [str sprintf("	{rank=same ")];
			for n = same	
				str = [str sprintf(" %i,", n)];
			end
			str(end) = sprintf(" ");
			str = [str sprintf("}\n")];
		end
		  
		str = [str sprintf("} \n\r")];
		fprintf(fid,str);
		fclose(fid);
		
		% Write xml file
		% Nodes
		
		fid=fopen( [fname '.xml'] , "w");
		str = [''];
		str = [str sprintf("<mc:Graph>\n")];
		str = [str sprintf(" <mc:NodeList>\n")];
		for i=1:nodes
			str = [str sprintf("	<mc:Node>\n")];
			str = [str sprintf("		<mc:id>%i</mc:id>\n",i)];
			str = [str sprintf("		<mc:name>P%i</mc:name>\n",i)];
			str = [str sprintf("		<mc:type>p</mc:type>\n")];
			str = [str sprintf("		<mc:rank>%i</mc:rank>\n",nranks(i))];
		end
	
		str = [str sprintf("	</mc:Node> \n\n")];
		str = [str sprintf(" </mc:NodeList>\n")];

		% Edges 
		str = [str sprintf(" <mc:EdgeList>\n")];
	
		for i=1:length(s) 
			str = [str sprintf("	<mc:Edge>\n")];
			str = [str sprintf("      	<mc:sourceId>%i</mc:sourceId>\n",s(i))];
			str = [str sprintf("       	<mc:targetId>%i</mc:targetId>\n",d(i))];
			str = [str sprintf("      	<mc:networkLoad>%i</mc:networkLoad>\n",loads(i))];
			str = [str sprintf("	</mc:Edge>\n")];
		end
			str = [str sprintf(" </mc:EdgeList>\n\n")];
		% Ranks
			str = [str sprintf(" <mc:RankList>\n")];

		for r=aranks
			str = [str sprintf("	<mc:RankGroup>\n")];
			rnodes = find(nranks == r);
			for n=rnodes
				str = [str sprintf("		<mc:Rank>%i</mc:Rank>\n", n)];
			end

			str = [str sprintf("	</mc:RankGroup>\n")];
		end 
			str = [str sprintf(" </mc:RankList>\n\n")];
			str = [str sprintf("</mc:Graph>\n")];

		fprintf(fid,str);
		fclose(fid);

		system(['dot -Teps ' fname '.gv -o ' fname '.eps']);
		system(['epspdf ' fname '.eps']);
%		system(['rm ' fname '.eps']); 

	end



