import { createClient } from '@supabase/supabase-js'

const supabaseUrl = import.meta.env.VITE_SUPABASE_URL
const supabaseAnonKey = import.meta.env.VITE_SUPABASE_ANON_KEY

if (!supabaseUrl || !supabaseAnonKey) {
  throw new Error('Missing Supabase environment variables')
}

export const supabase = createClient(supabaseUrl, supabaseAnonKey)

export type Database = {
  public: {
    Tables: {
      users: {
        Row: {
          id: string
          email: string
          name: string
          avatar_url: string | null
          created_at: string | null
        }
        Insert: {
          id?: string
          email: string
          name?: string
          avatar_url?: string | null
          created_at?: string | null
        }
        Update: {
          id?: string
          email?: string
          name?: string
          avatar_url?: string | null
          created_at?: string | null
        }
      }
      labs: {
        Row: {
          id: string
          name: string
          description: string | null
          owner_id: string
          created_at: string | null
        }
        Insert: {
          id?: string
          name?: string
          description?: string | null
          owner_id: string
          created_at?: string | null
        }
        Update: {
          id?: string
          name?: string
          description?: string | null
          owner_id?: string
          created_at?: string | null
        }
      }
      lab_members: {
        Row: {
          id: string
          user_id: string
          lab_id: string
          role: 'admin' | 'member'
          joined_at: string | null
          member_name: string
        }
        Insert: {
          id?: string
          user_id: string
          lab_id: string
          role?: 'admin' | 'member'
          joined_at?: string | null
          member_name?: string
        }
        Update: {
          id?: string
          user_id?: string
          lab_id?: string
          role?: 'admin' | 'member'
          joined_at?: string | null
          member_name?: string
        }
      }
      invitations: {
        Row: {
          id: string
          lab_id: string
          invited_email: string
          invited_by: string
          status: 'pending' | 'accepted' | 'rejected'
          created_at: string | null
          expires_at: string | null
          lab_name: string
          invited_by_name: string
        }
        Insert: {
          id?: string
          lab_id: string
          invited_email: string
          invited_by: string
          status?: 'pending' | 'accepted' | 'rejected'
          created_at?: string | null
          expires_at?: string | null
          lab_name?: string
          invited_by_name?: string
        }
        Update: {
          id?: string
          lab_id?: string
          invited_email?: string
          invited_by?: string
          status?: 'pending' | 'accepted' | 'rejected'
          created_at?: string | null
          expires_at?: string | null
          lab_name?: string
          invited_by_name?: string
        }
      }
      task_lists: {
        Row: {
          id: string
          lab_id: string
          name: string
          description: string | null
          color: string
          created_by: string
          created_at: string | null
          updated_at: string | null
        }
        Insert: {
          id?: string
          lab_id: string
          name?: string
          description?: string | null
          color?: string
          created_by: string
          created_at?: string | null
          updated_at?: string | null
        }
        Update: {
          id?: string
          lab_id?: string
          name?: string
          description?: string | null
          color?: string
          created_by?: string
          created_at?: string | null
          updated_at?: string | null
        }
      }
      tasks: {
        Row: {
          id: string
          task_list_id: string
          title: string
          description: string | null
          completed: boolean | null
          completed_by: string | null
          completed_at: string | null
          on_it_by: string | null
          on_it_at: string | null
          assigned_to: string | null
          created_by: string
          created_at: string | null
          updated_at: string | null
          priority: string | null
          due_date: string | null
        }
        Insert: {
          id?: string
          task_list_id: string
          title?: string
          description?: string | null
          completed?: boolean | null
          completed_by?: string | null
          completed_at?: string | null
          on_it_by?: string | null
          on_it_at?: string | null
          assigned_to?: string | null
          created_by: string
          created_at?: string | null
          updated_at?: string | null
          priority?: string | null
          due_date?: string | null
        }
        Update: {
          id?: string
          task_list_id?: string
          title?: string
          description?: string | null
          completed?: boolean | null
          completed_by?: string | null
          completed_at?: string | null
          on_it_by?: string | null
          on_it_at?: string | null
          assigned_to?: string | null
          created_by?: string
          created_at?: string | null
          updated_at?: string | null
          priority?: string | null
          due_date?: string | null
        }
      }
      task_activities: {
        Row: {
          id: string
          task_id: string
          lab_id: string
          user_id: string
          activity_type: string
          message: string
          created_at: string | null
          read_by: any | null
        }
        Insert: {
          id?: string
          task_id: string
          lab_id: string
          user_id: string
          activity_type: string
          message: string
          created_at?: string | null
          read_by?: any | null
        }
        Update: {
          id?: string
          task_id?: string
          lab_id?: string
          user_id?: string
          activity_type?: string
          message?: string
          created_at?: string | null
          read_by?: any | null
        }
      }
    }
  }
}