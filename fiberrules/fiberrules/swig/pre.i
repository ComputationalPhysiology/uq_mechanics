/* -*- C -*- */
// Copyright (C) 2012 Johan Hake
//
// This file is part of FIBERRULES.
//
// FIBERRULES is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// FIBERRULES is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with FIBERRULES. If not, see <http://www.gnu.org/licenses/>.

%ignore Element::indices;
%rename (_readElemFile) readElemFile;

// Reduce vector interface
%ignore std::vector<Element>::capacity;
%ignore std::vector<Element>::assign;
%ignore std::vector<Element>::erase;
%ignore std::vector<Element>::iterator;
%ignore std::vector<Element>::empty;
%ignore std::vector<Element>::begin;
%ignore std::vector<Element>::end;
%ignore std::vector<Element>::rbegin;
%ignore std::vector<Element>::rend;
%ignore std::vector<Element>::front;
%ignore std::vector<Element>::get_allocator;
%ignore std::vector<Element>::pop;
%ignore std::vector<Element>::pop_back;
%ignore std::vector<Element>::swap;
%ignore std::vector<Element>::insert;
%ignore std::vector<Element>::push_back;

