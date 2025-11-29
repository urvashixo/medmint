/*
  # Remove All RLS Policies

  This migration removes all Row Level Security policies from all tables.
  Security will be managed at the application level instead.

  1. Disable RLS on all tables
  2. Drop all existing policies
  3. Remove helper functions

  Note: This approach moves security responsibility to the application layer,
  which can be more flexible and easier to debug than complex RLS policies.
*/

-- Disable Row Level Security on all tables
ALTER TABLE users DISABLE ROW LEVEL SECURITY;
ALTER TABLE labs DISABLE ROW LEVEL SECURITY;
ALTER TABLE lab_members DISABLE ROW LEVEL SECURITY;
ALTER TABLE whiteboards DISABLE ROW LEVEL SECURITY;
ALTER TABLE reports DISABLE ROW LEVEL SECURITY;
ALTER TABLE compounds DISABLE ROW LEVEL SECURITY;
ALTER TABLE proteins DISABLE ROW LEVEL SECURITY;
ALTER TABLE chat_messages DISABLE ROW LEVEL SECURITY;
ALTER TABLE todos DISABLE ROW LEVEL SECURITY;

-- Drop all existing policies on users
DROP POLICY IF EXISTS "Users can read own profile" ON users;
DROP POLICY IF EXISTS "Users can update own profile" ON users;
DROP POLICY IF EXISTS "Users can insert own profile" ON users;

-- Drop all existing policies on labs
DROP POLICY IF EXISTS "Lab owners can manage their labs" ON labs;
DROP POLICY IF EXISTS "Lab members can read labs where they are members" ON labs;
DROP POLICY IF EXISTS "Users can read labs where they are members" ON labs;

-- Drop all existing policies on lab_members
DROP POLICY IF EXISTS "Lab admins can manage lab members" ON lab_members;
DROP POLICY IF EXISTS "Users can read lab memberships" ON lab_members;
DROP POLICY IF EXISTS "Lab owners can add members" ON lab_members;
DROP POLICY IF EXISTS "Lab owners can read all lab memberships" ON lab_members;
DROP POLICY IF EXISTS "Lab owners can remove members" ON lab_members;
DROP POLICY IF EXISTS "Lab owners can update member roles" ON lab_members;
DROP POLICY IF EXISTS "Users can leave labs" ON lab_members;
DROP POLICY IF EXISTS "Users can read own lab memberships" ON lab_members;
DROP POLICY IF EXISTS "Lab owners can manage members" ON lab_members;
DROP POLICY IF EXISTS "Lab admins can manage members" ON lab_members;

-- Drop all existing policies on whiteboards
DROP POLICY IF EXISTS "Lab members can manage whiteboards" ON whiteboards;

-- Drop all existing policies on reports
DROP POLICY IF EXISTS "Lab members can read reports" ON reports;
DROP POLICY IF EXISTS "Lab members can create reports" ON reports;
DROP POLICY IF EXISTS "Report creators can update their reports" ON reports;

-- Drop all existing policies on compounds
DROP POLICY IF EXISTS "Lab members can read compounds" ON compounds;
DROP POLICY IF EXISTS "Lab members can create compounds" ON compounds;
DROP POLICY IF EXISTS "Compound creators can update their compounds" ON compounds;

-- Drop all existing policies on proteins
DROP POLICY IF EXISTS "Lab members can read proteins" ON proteins;
DROP POLICY IF EXISTS "Lab members can create proteins" ON proteins;
DROP POLICY IF EXISTS "Protein creators can update their proteins" ON proteins;

-- Drop all existing policies on chat_messages
DROP POLICY IF EXISTS "Lab members can read chat messages" ON chat_messages;
DROP POLICY IF EXISTS "Lab members can send chat messages" ON chat_messages;

-- Drop all existing policies on todos
DROP POLICY IF EXISTS "Lab members can read todos" ON todos;
DROP POLICY IF EXISTS "Lab members can create todos" ON todos;
DROP POLICY IF EXISTS "Todo creators can update their todos" ON todos;
DROP POLICY IF EXISTS "Lab admins can update any todos" ON todos;

-- Drop helper functions that are no longer needed
DROP FUNCTION IF EXISTS is_lab_owner(uuid, uuid);

-- Note: The triggers and trigger functions for lab owner auto-add and whiteboard timestamps
-- are kept as they provide useful functionality without security implications